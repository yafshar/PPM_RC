      !-------------------------------------------------------------------------
      !  Module       :               ppm_rc_module_topologicalnumber
      !-------------------------------------------------------------------------
      !
      !  Purpose      : provides some code to compute the topological numbers
      !                 introduced by Giles Bertrand and Gregoire Malandin.
      !
      !             COMPATIBLE CONNECTIVITIES / TOPOLOGICAL NUMBERS
      !
      !             Regular connectivity       Cell Connectivity
      !      1:     4, in two dimensions            (2,1)
      !      2:     8, in two dimensions            (2,0)
      !      3:     6, in three dimensions          (3,2)
      !      3:    18, in three dimensions          (3,1)
      !      3:    26, in three dimensions          (3,0)
      !      1:            (6+,18)
      !      2:            (18,6+)
      !      3:            (6,26)
      !      4:            (26,6)
      !
      !  Remarks      :
      !
      !  References   : G. Bertrand and G. Malandain : "A new characterization
      !                 of three-dimensional simple points",
      !                 Pattern Recognition Letters; 15:169--175; 1994.
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      MODULE ppm_rc_module_topologicalnumber
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_rc_module_global
        IMPLICIT NONE

        PRIVATE

        !!! Direction arrays for jumping around in mesh space
        !!! based on topology connectivity
        INTEGER, DIMENSION(2,4),  PARAMETER :: Offset_4 =  &
        & RESHAPE((/       0,-1,                           &
        &           -1, 0,       1, 0,                     &
        &                  0, 1      /),           (/2, 4/))
        !!! 4-neighborhood in two-dimension
        INTEGER, DIMENSION(2,8),  PARAMETER :: Offset_2d = &
        & RESHAPE((/-1,-1, 0,-1, 1,-1,                     &
        &           -1, 0,       1, 0,                     &
        &           -1, 1, 0, 1, 1, 1/),           (/2, 8/))
        !!! 8-neighborhood in two-dimension
        INTEGER, DIMENSION(3,6),  PARAMETER :: Offset_6 =  &
        & RESHAPE((/           0, 0,-1,                    &
        &                      0,-1, 0,                    &
        &           -1, 0, 0,           1, 0, 0,           &
        &                      0, 1, 0,                    &
        &                      0, 0, 1         /), (/3, 6/))
        !!! 6-neighborhood in three-dimension
        INTEGER, DIMENSION(3,18), PARAMETER :: Offset_18 = &
        & RESHAPE((/           0,-1,-1,                    &
        &           -1, 0,-1,  0, 0,-1, 1, 0,-1,           &
        &                      0, 1,-1,                    &
        &           -1,-1, 0,  0,-1, 0, 1,-1, 0,           &
        &           -1, 0, 0,           1, 0, 0,           &
        &           -1, 1, 0,  0, 1, 0, 1, 1, 0,           &
        &                      0,-1, 1,                    &
        &           -1, 0, 1,  0, 0, 1, 1, 0, 1,           &
        &                      0, 1, 1         /), (/3,18/))
        !!! 18-neighborhood in three-dimension
        INTEGER, DIMENSION(3,26), PARAMETER :: Offset_3d = &
        & RESHAPE((/-1,-1,-1,  0,-1,-1, 1,-1,-1,           &
        &           -1, 0,-1,  0, 0,-1, 1, 0,-1,           &
        &           -1, 1,-1,  0, 1,-1, 1, 1,-1,           &
        &           -1,-1, 0,  0,-1, 0, 1,-1, 0,           &
        &           -1, 0, 0,           1, 0, 0,           &
        &           -1, 1, 0,  0, 1, 0, 1, 1, 0,           &
        &           -1,-1, 1,  0,-1, 1, 1,-1, 1,           &
        &           -1, 0, 1,  0, 0, 1, 1, 0, 1,           &
        &           -1, 1, 1,  0, 1, 1, 1, 1, 1/), (/3,26/))
        !!! 26-neighborhood in three-dimension

!         PUBLIC :: Offset_4,Offset_2d,Offset_6,Offset_18,Offset_3d

        !----------------------------------------------------------------------
        !
        !----------------------------------------------------------------------
#include "./topological_number/ppm_rc_toponumtypedef.f"

        TYPE(TopologicalNumberImageFunction) :: TopologicalNumberFunction

        TYPE(Connectivity), POINTER :: FG_ConnectivityType => NULL()
        TYPE(Connectivity), POINTER :: BG_ConnectivityType => NULL()

        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE IsSingleFGPoint
          MODULE PROCEDURE IsSingleFGPoint_2d
          MODULE PROCEDURE IsSingleFGPoint__2d
          MODULE PROCEDURE IsSingleFGPoint_3d
          MODULE PROCEDURE IsSingleFGPoint__3d
        END INTERFACE

        INTERFACE IsEnclosedByLabel_FGConnectivity
          MODULE PROCEDURE IsEnclosedByLabel_FGConnectivity_2d
          MODULE PROCEDURE IsEnclosedByLabel_FGConnectivity__2d
          MODULE PROCEDURE IsEnclosedByLabel_FGConnectivity_3d
          MODULE PROCEDURE IsEnclosedByLabel_FGConnectivity__3d
        END INTERFACE

        INTERFACE IsEnclosedByLabel_BGConnectivity
          MODULE PROCEDURE IsEnclosedByLabel_BGConnectivity_2d
          MODULE PROCEDURE IsEnclosedByLabel_BGConnectivity__2d
          MODULE PROCEDURE IsEnclosedByLabel_BGConnectivity_3d
          MODULE PROCEDURE IsEnclosedByLabel_BGConnectivity__3d
        END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: FGandBGTopoNbPairType
        PUBLIC :: TopologicalNumberFunction
        PUBLIC :: FG_ConnectivityType
        PUBLIC :: BG_ConnectivityType
        PUBLIC :: Connectivity
        PUBLIC :: FGNeighborhoodConnectivity
        PUBLIC :: BGNeighborhoodConnectivity
        PUBLIC :: BackgroundConnectivity

        PUBLIC :: IsSingleFGPoint
        PUBLIC :: IsEnclosedByLabel_FGConnectivity
        PUBLIC :: IsEnclosedByLabel_BGConnectivity

      CONTAINS

#include "./topological_number/ppm_rc_toponumtypeproc.f"

      END MODULE ppm_rc_module_topologicalnumber
