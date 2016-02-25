        ! Connectivity information.
        !
        ! This TYPE describes the k-neighbors of a point in a digital space of
        ! dimension 2 or 3.
        ! The definition used for this connectivity is the one using cell
        ! decomposition, as defined by Malandain in "On Topology in Multidimensional
        ! Discrete Spaces", ftp://ftp.inria.fr/INRIA/publication/RR/RR-2098.ps.gz
        !
        ! The advantage of this definition is that it is consistent and
        ! extensible in n-dimensions.
        !
        ! In 2D, 4- and 8-connectivity are respectively corresponding to 1- and
        ! 0-connectivities. In 3D, 6-, 18- and 26-connectivity are respectively
        ! corresponding to 2-, 1- and 0-connectivity.
        TYPE :: Connectivity
          !!! Connectivity
          INTEGER                              :: ID=3
          !!! Connectivity type
          INTEGER                              :: Dimension=3
          !!! The dimension of the space
          INTEGER                              :: CellDimension=2
          !!! The dimension of the cell based on cell decomposition
          !!! which will define the type of neighborhood
          INTEGER                              :: NeighborhoodSize=26
          !!! Full number of neighbors independent on the connectivity type
          !!! (only is dependent on the dimension of space)
          INTEGER                              :: NumberOfNeighbors=6
          !!! Number of neighbors with respect to the connectivity type
          INTEGER, DIMENSION(:),   ALLOCATABLE :: NeighborsOffsets
          !!! Neighbors as offsets
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: NeighborsPoints

          TYPE(Connectivity),      POINTER     :: FGNeighborhood => NULL()
          TYPE(Connectivity),      POINTER     :: BGNeighborhood => NULL()
          TYPE(Connectivity),      POINTER     :: Background     => NULL()

        CONTAINS
          ! Constructor for connectivity structure
          PROCEDURE :: create => connectivity_create
          PROCEDURE :: destroy => connectivity_destroy

          ! returns the full number of neighbors independent
          ! on the connectivity type (only dependent on the dimension)
          PROCEDURE :: GetNeighborhoodSize
          ! returns the number of neighbors with respect to the
          ! connectivity type.
          PROCEDURE :: GetNumberOfNeighbors

          PROCEDURE :: GetNeighborsOffsets
!           PROCEDURE :: GetNeighborsPoints

          ! Test if two points are neighbors
          PROCEDURE :: PointsAreNeighbors
          !TODO
          !I should implement the procedure
!           PROCEDURE :: OffsetsAreNeighbors
          GENERIC   :: AreNeighbors =>     &
          &            PointsAreNeighbors  !, &
!           &            OffsetsAreNeighbors

          !Test if a point is a neighbor of 0
          PROCEDURE :: PointIsInNeighborhood
          PROCEDURE :: OffsetIsInNeighborhood
          GENERIC   :: IsInNeighborhood =>    &
          &            PointIsInNeighborhood, &
          &            OffsetIsInNeighborhood

          !Convert an offset to a point, in a 2x2 or 3x3x3 cube
          PROCEDURE :: OffsetToPoint
          !Convert a point to an offset, in a 2x2 or 3x3x3 cube
          PROCEDURE :: PointToOffset

        END TYPE Connectivity

        ! Helper to determine the neighborhood connectivity to use when
        ! computing topological numbers.
        ! Connectivity : the image connectivity
        INTERFACE FGNeighborhoodConnectivity
          MODULE PROCEDURE FGNeighborhoodConnectivity
        END INTERFACE FGNeighborhoodConnectivity

        INTERFACE BGNeighborhoodConnectivity
          MODULE PROCEDURE BGNeighborhoodConnectivity
        END INTERFACE BGNeighborhoodConnectivity

        ! Helper to determine the background connectivity
        ! FGConnectivity : the image connectivity
        INTERFACE BackgroundConnectivity
          MODULE PROCEDURE BackgroundConnectivity
        END INTERFACE BackgroundConnectivity

        !  For each point in the n'-neighborhood of 0, characterizes which of
        !  its surrounding points are n-neighbors and belong to the [-1, 1] cube.
        !
        !  Connectivity : the image connectivity, (n in the description)
        !  NeighborhoodConnectivity : the neighborhood connectivity, (n' in the description)
        !
        !  In 3D, use (6,18) or (26,26) for (n, n') to get significant result for
        !  topological numbers.
        !
        !  The criterion returns true for (p1, p2) if p1 is n'-adjacent to 0 and
        !  (p1, p1+p2) are n-adjacent and (p1+p2) is in [-1, 1]^3.
        !!! constructor
        INTERFACE UnitCubeNeighbors
          MODULE PROCEDURE constructor
        END INTERFACE UnitCubeNeighbors

        ! Counting the number of connected components restricted in a
        ! unit cube. This is used for topological number computation.
        TYPE UnitCubeCCCounter
          TYPE(Connectivity),      POINTER     :: TConnectivity => NULL()
          TYPE(Connectivity),      POINTER     :: TNeighborhoodConnectivity => NULL()

          LOGICAL, DIMENSION(:),   ALLOCATABLE :: m_Image
          LOGICAL, DIMENSION(:),   ALLOCATABLE :: ConnectivityTest
          LOGICAL, DIMENSION(:,:), POINTER     :: UnitCubeNeighbors => NULL()

        CONTAINS
          ! create structure
          PROCEDURE :: create => UnitCubeCCCounter_create
          ! destroy structure
          PROCEDURE :: destroy => UnitCubeCCCounter_destroy

          PROCEDURE :: SetImage

          PROCEDURE :: ConnectedComponentCounter

        END TYPE UnitCubeCCCounter

        TYPE FGandBGTopoNbPairType
          INTEGER :: FGNumber
          INTEGER :: BGNumber
        END TYPE FGandBGTopoNbPairType

        TYPE ForegroundTopologicalNumberType
          INTEGER                     :: label
          TYPE(FGandBGTopoNbPairType) :: FGandBGTopoNbPair
        CONTAINS
          PROCEDURE :: ForegroundTopologicalNumberType_create1
          PROCEDURE :: ForegroundTopologicalNumberType_create2
          GENERIC   :: create =>                                &
          &            ForegroundTopologicalNumberType_create1, &
          &            ForegroundTopologicalNumberType_create2
        END TYPE ForegroundTopologicalNumberType

        ! Compute the topological numbers of an image at a given index.
        !
        ! Topological numbers characterize the topological properties of a point. They
        ! are defined in an article by G. Bertrand and G. Malandain : "A new
        ! characterization of three-dimensional simple points"; Pattern Recognition
        ! Letters; 15:169--175; 1994.
        TYPE TopologicalNumberImageFunction
          TYPE(Connectivity),      POINTER     :: TFGConnectivity => NULL()
          TYPE(Connectivity),      POINTER     :: TBGConnectivity => NULL()

          TYPE(UnitCubeCCCounter)              :: ForegroundUnitCubeCCCounter
          TYPE(UnitCubeCCCounter)              :: BackgroundUnitCubeCCCounter

          INTEGER, DIMENSION(:),   ALLOCATABLE :: DataSubImage
          ! cached input image

          LOGICAL                              :: ComputeForegroundTN=.TRUE.
          LOGICAL                              :: ComputeBackgroundTN=.TRUE.
          ! These two members allow to selectively compute the topological
          ! numbers for the background and the foreground. They are both set to true
          ! during the construction of the object.

          LOGICAL, DIMENSION(:),   ALLOCATABLE :: SubImage
          !binary subimage (switching fg/bg)

        CONTAINS

          ! create structure
          PROCEDURE :: create => TopologicalNumberImageFunction_create
          ! destroy structure
          PROCEDURE :: destroy => TopologicalNumberImageFunction_destroy

          PROCEDURE :: readDataSubImagePoint
          PROCEDURE :: readDataSubImage_2d
          PROCEDURE :: readDataSubImage_3d
          PROCEDURE :: readDataSubImage__2d
          PROCEDURE :: readDataSubImage__3d
          GENERIC   :: readDataSubImage =>    &
          &            readDataSubImagePoint, &
          &            readDataSubImage_2d,   &
          &            readDataSubImage_3d,   &
          &            readDataSubImage__2d,  &
          &            readDataSubImage__3d

          PROCEDURE :: EvaluateAtIndex1_2d
          PROCEDURE :: EvaluateAtIndex1_3d
          PROCEDURE :: EvaluateAtIndex1__2d
          PROCEDURE :: EvaluateAtIndex1__3d
          PROCEDURE :: EvaluateAtIndex2
          PROCEDURE :: EvaluateAtIndex3
          GENERIC   :: EvaluateAtIndex =>    &
          &            EvaluateAtIndex1_2d,  &
          &            EvaluateAtIndex1_3d,  &
          &            EvaluateAtIndex1__2d, &
          &            EvaluateAtIndex1__3d, &
          &            EvaluateAtIndex2,     &
          &            EvaluateAtIndex3

          PROCEDURE :: EvaluateFGTNOfLabelAtIndex1_2d
          PROCEDURE :: EvaluateFGTNOfLabelAtIndex1_3d
          PROCEDURE :: EvaluateFGTNOfLabelAtIndex1__2d
          PROCEDURE :: EvaluateFGTNOfLabelAtIndex1__3d
          PROCEDURE :: EvaluateFGTNOfLabelAtIndex2
          PROCEDURE :: EvaluateFGTNOfLabelAtIndex3
          GENERIC   :: EvaluateFGTNOfLabelAtIndex =>    &
          &            EvaluateFGTNOfLabelAtIndex1_2d,  &
          &            EvaluateFGTNOfLabelAtIndex1_3d,  &
          &            EvaluateFGTNOfLabelAtIndex1__2d, &
          &            EvaluateFGTNOfLabelAtIndex1__3d, &
          &            EvaluateFGTNOfLabelAtIndex2,     &
          &            EvaluateFGTNOfLabelAtIndex3

          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex1_2d
          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex1_3d
          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex1__2d
          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex1__3d
          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex2
          PROCEDURE :: EvaluateAdjacentRegionsFGTNAtIndex3
          GENERIC   :: EvaluateAdjacentRegionsFGTNAtIndex =>    &
          &            EvaluateAdjacentRegionsFGTNAtIndex1_2d,  &
          &            EvaluateAdjacentRegionsFGTNAtIndex1_3d,  &
          &            EvaluateAdjacentRegionsFGTNAtIndex1__2d, &
          &            EvaluateAdjacentRegionsFGTNAtIndex1__3d, &
          &            EvaluateAdjacentRegionsFGTNAtIndex2,     &
          &            EvaluateAdjacentRegionsFGTNAtIndex3

        END TYPE TopologicalNumberImageFunction
