        ! Constructor for connectivity structure
        SUBROUTINE Connectivity_create(this,dim,connectivity_num,info)

          IMPLICIT NONE

          CLASS(Connectivity)    :: this

          INTEGER, INTENT(IN   ) :: dim
          INTEGER, INTENT(IN   ) :: connectivity_num
          INTEGER, INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("Connectivity_create")

          SELECT CASE (connectivity_num)
          CASE (4)
             IF (dim.NE.2) THEN
                fail("connectivity has conflict with the dimension")
             ENDIF
             this%ID=1
             this%Dimension=2
             this%CellDimension=1
             ALLOCATE(this%NeighborsOffsets(4),STAT=info)
             or_fail_alloc("this%NeighborsOffsets")
             ALLOCATE(this%NeighborsPoints(2,4),SOURCE=Offset_4,STAT=info)
             or_fail_alloc("this%NeighborsPoints")

          CASE (8)
             IF (dim.NE.2) THEN
                fail("connectivity has conflict with the dimension")
             ENDIF
             this%ID=2
             this%Dimension=2
             this%CellDimension=0
             ALLOCATE(this%NeighborsOffsets(8),STAT=info)
             or_fail_alloc("this%NeighborsOffsets")
             ALLOCATE(this%NeighborsPoints(2,8),SOURCE=Offset_2d,STAT=info)
             or_fail_alloc("this%NeighborsPoints")

          CASE (6)
             IF (dim.NE.3) THEN
                fail("connectivity has conflict with the dimension")
             ENDIF
             this%ID=3
             this%Dimension=3
             this%CellDimension=2
             ALLOCATE(this%NeighborsOffsets(6),STAT=info)
             or_fail_alloc("this%NeighborsOffsets")
             ALLOCATE(this%NeighborsPoints(3,6),SOURCE=Offset_6,STAT=info)
             or_fail_alloc("this%NeighborsPoints")

          CASE (18)
             IF (dim.NE.3) THEN
                fail("connectivity has conflict with the dimension")
             ENDIF
             this%ID=4
             this%Dimension=3
             this%CellDimension=1
             ALLOCATE(this%NeighborsOffsets(18),STAT=info)
             or_fail_alloc("this%NeighborsOffsets")
             ALLOCATE(this%NeighborsPoints(3,18),SOURCE=Offset_18,STAT=info)
             or_fail_alloc("this%NeighborsPoints")

          CASE (26)
             IF (dim.NE.3) THEN
                fail("connectivity has conflict with the dimension")
             ENDIF
             this%ID=5
             this%Dimension=3
             this%CellDimension=0
             ALLOCATE(this%NeighborsOffsets(26),STAT=info)
             or_fail_alloc("this%NeighborsOffsets")
             ALLOCATE(this%NeighborsPoints(3,26),SOURCE=Offset_3d,STAT=info)
             or_fail_alloc("this%NeighborsPoints")

          CASE DEFAULT
             fail("dimension and connectivity does not match")

          END SELECT

          CALL this%GetNeighborhoodSize()
          CALL this%GetNumberOfNeighbors()
          CALL this%GetNeighborsOffsets()

          end_subroutine()

        END SUBROUTINE Connectivity_create

        ! Constructor for connectivity structure
        SUBROUTINE connectivity_destroy(this,info)

          IMPLICIT NONE

          CLASS(Connectivity)    :: this

          INTEGER, INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          start_subroutine("connectivity_destroy")

          DEALLOCATE(this%NeighborsOffsets,STAT=info)
          or_fail_dealloc("this%NeighborsOffsets")

          DEALLOCATE(this%NeighborsPoints,STAT=info)
          or_fail_dealloc("this%NeighborsPoints")

          IF (ASSOCIATED(this%FGNeighborhood)) THEN
             ASSOCIATE (FGN => this%FGNeighborhood)
                CALL FGN%destroy(info)
                or_fail("Failed to destroy this%FGNeighborhood!")
             END ASSOCIATE
             DEALLOCATE(this%FGNeighborhood,STAT=info)
             or_fail("Failed to deallocate FGNeighborhood!")
          ENDIF

          IF (ASSOCIATED(this%BGNeighborhood)) THEN
             ASSOCIATE (BGN => this%BGNeighborhood)
                CALL BGN%destroy(info)
                or_fail("Failed to destroy BGNeighborhood!")
             END ASSOCIATE
             DEALLOCATE(this%BGNeighborhood,STAT=info)
             or_fail("Failed to deallocate BGNeighborhood!")
          ENDIF

          IF (ASSOCIATED(this%Background)) THEN
             ASSOCIATE (BGC => this%Background)
                CALL BGC%destroy(info)
                or_fail("Failed to destroy this%Background!")
             END ASSOCIATE
             DEALLOCATE(this%Background,STAT=info)
             or_fail("Failed to deallocate Background!")
          ENDIF

          NULLIFY(this%FGNeighborhood)
          NULLIFY(this%BGNeighborhood)
          NULLIFY(this%Background)

          end_subroutine()

        END SUBROUTINE connectivity_destroy

        ! returns the full number of neighbors independent
        ! on the connectivity type (only dependent on the dimension)
        SUBROUTINE GetNeighborhoodSize(this)

          IMPLICIT NONE

          CLASS(Connectivity) :: this

          SELECT CASE (this%ID)
          CASE (1,2)
             this%NeighborhoodSize=8

          CASE DEFAULT
             this%NeighborhoodSize=26

          END SELECT
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE GetNeighborhoodSize

        ! returns the number of neighbors with respect to the
        ! connectivity type.
        SUBROUTINE GetNumberOfNeighbors(this)
          IMPLICIT NONE

          CLASS(Connectivity) :: this

          SELECT CASE (this%ID)
          !4, in two dimensions
          CASE (1)
             this%NumberOfNeighbors=4
          !8, in two dimensions
          CASE (2)
             this%NumberOfNeighbors=8
          !6, in three dimensions
          CASE (3)
             this%NumberOfNeighbors=6
          !18, in three dimensions
          CASE (4)
             this%NumberOfNeighbors=18
          !26, in three dimensions
          CASE (5)
             this%NumberOfNeighbors=26

          END SELECT
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE GetNumberOfNeighbors

        SUBROUTINE GetNeighborsOffsets(this)

          IMPLICIT NONE

          CLASS(Connectivity)    :: this

          SELECT CASE (this%ID)
          !4, in two dimensions
          CASE (1)
             this%NeighborsOffsets(1)=this%PointToOffset(Offset_4(:,1))
             this%NeighborsOffsets(2)=this%PointToOffset(Offset_4(:,2))
             this%NeighborsOffsets(3)=this%PointToOffset(Offset_4(:,3))
             this%NeighborsOffsets(4)=this%PointToOffset(Offset_4(:,4))

          !8, in two dimensions
          CASE (2)
             this%NeighborsOffsets(1)=this%PointToOffset(Offset_2d(:,1))
             this%NeighborsOffsets(2)=this%PointToOffset(Offset_2d(:,2))
             this%NeighborsOffsets(3)=this%PointToOffset(Offset_2d(:,3))
             this%NeighborsOffsets(4)=this%PointToOffset(Offset_2d(:,4))
             this%NeighborsOffsets(5)=this%PointToOffset(Offset_2d(:,5))
             this%NeighborsOffsets(6)=this%PointToOffset(Offset_2d(:,6))
             this%NeighborsOffsets(7)=this%PointToOffset(Offset_2d(:,7))
             this%NeighborsOffsets(8)=this%PointToOffset(Offset_2d(:,8))

          !6, in three dimensions
          CASE (3)
             this%NeighborsOffsets(1)=this%PointToOffset(Offset_6(:,1))
             this%NeighborsOffsets(2)=this%PointToOffset(Offset_6(:,2))
             this%NeighborsOffsets(3)=this%PointToOffset(Offset_6(:,3))
             this%NeighborsOffsets(4)=this%PointToOffset(Offset_6(:,4))
             this%NeighborsOffsets(5)=this%PointToOffset(Offset_6(:,5))
             this%NeighborsOffsets(6)=this%PointToOffset(Offset_6(:,6))

          !18, in three dimensions
          CASE (4)
             this%NeighborsOffsets( 1)=this%PointToOffset(Offset_18(:, 1))
             this%NeighborsOffsets( 2)=this%PointToOffset(Offset_18(:, 2))
             this%NeighborsOffsets( 3)=this%PointToOffset(Offset_18(:, 3))
             this%NeighborsOffsets( 4)=this%PointToOffset(Offset_18(:, 4))
             this%NeighborsOffsets( 5)=this%PointToOffset(Offset_18(:, 5))
             this%NeighborsOffsets( 6)=this%PointToOffset(Offset_18(:, 6))
             this%NeighborsOffsets( 7)=this%PointToOffset(Offset_18(:, 7))
             this%NeighborsOffsets( 8)=this%PointToOffset(Offset_18(:, 8))
             this%NeighborsOffsets( 9)=this%PointToOffset(Offset_18(:, 9))
             this%NeighborsOffsets(10)=this%PointToOffset(Offset_18(:,10))
             this%NeighborsOffsets(11)=this%PointToOffset(Offset_18(:,11))
             this%NeighborsOffsets(12)=this%PointToOffset(Offset_18(:,12))
             this%NeighborsOffsets(13)=this%PointToOffset(Offset_18(:,13))
             this%NeighborsOffsets(14)=this%PointToOffset(Offset_18(:,14))
             this%NeighborsOffsets(15)=this%PointToOffset(Offset_18(:,15))
             this%NeighborsOffsets(16)=this%PointToOffset(Offset_18(:,16))
             this%NeighborsOffsets(17)=this%PointToOffset(Offset_18(:,17))
             this%NeighborsOffsets(18)=this%PointToOffset(Offset_18(:,18))

          !26, in three dimensions
          CASE (5)
             this%NeighborsOffsets( 1)=this%PointToOffset(Offset_3d(:, 1))
             this%NeighborsOffsets( 2)=this%PointToOffset(Offset_3d(:, 2))
             this%NeighborsOffsets( 3)=this%PointToOffset(Offset_3d(:, 3))
             this%NeighborsOffsets( 4)=this%PointToOffset(Offset_3d(:, 4))
             this%NeighborsOffsets( 5)=this%PointToOffset(Offset_3d(:, 5))
             this%NeighborsOffsets( 6)=this%PointToOffset(Offset_3d(:, 6))
             this%NeighborsOffsets( 7)=this%PointToOffset(Offset_3d(:, 7))
             this%NeighborsOffsets( 8)=this%PointToOffset(Offset_3d(:, 8))
             this%NeighborsOffsets( 9)=this%PointToOffset(Offset_3d(:, 9))
             this%NeighborsOffsets(10)=this%PointToOffset(Offset_3d(:,10))
             this%NeighborsOffsets(11)=this%PointToOffset(Offset_3d(:,11))
             this%NeighborsOffsets(12)=this%PointToOffset(Offset_3d(:,12))
             this%NeighborsOffsets(13)=this%PointToOffset(Offset_3d(:,13))
             this%NeighborsOffsets(14)=this%PointToOffset(Offset_3d(:,14))
             this%NeighborsOffsets(15)=this%PointToOffset(Offset_3d(:,15))
             this%NeighborsOffsets(16)=this%PointToOffset(Offset_3d(:,16))
             this%NeighborsOffsets(17)=this%PointToOffset(Offset_3d(:,17))
             this%NeighborsOffsets(18)=this%PointToOffset(Offset_3d(:,18))
             this%NeighborsOffsets(19)=this%PointToOffset(Offset_3d(:,19))
             this%NeighborsOffsets(20)=this%PointToOffset(Offset_3d(:,20))
             this%NeighborsOffsets(21)=this%PointToOffset(Offset_3d(:,21))
             this%NeighborsOffsets(22)=this%PointToOffset(Offset_3d(:,22))
             this%NeighborsOffsets(23)=this%PointToOffset(Offset_3d(:,23))
             this%NeighborsOffsets(24)=this%PointToOffset(Offset_3d(:,24))
             this%NeighborsOffsets(25)=this%PointToOffset(Offset_3d(:,25))
             this%NeighborsOffsets(26)=this%PointToOffset(Offset_3d(:,26))

          END SELECT
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE GetNeighborsOffsets

        LOGICAL FUNCTION PointsAreNeighbors(this,coord1,coord2)

          IMPLICIT NONE

          CLASS(Connectivity)                  :: this

          INTEGER, DIMENSION(:), INTENT(IN   ) :: coord1
          INTEGER, DIMENSION(:), INTENT(IN   ) :: coord2

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(SIZE(coord1)) :: difference

          difference=coord1-coord2

          PointsAreNeighbors=this%IsInNeighborhood(difference)

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION PointsAreNeighbors

        LOGICAL FUNCTION PointIsInNeighborhood(this,coord)

          IMPLICIT NONE

          CLASS(Connectivity)                  :: this

          INTEGER, DIMENSION(:), INTENT(IN   ) :: coord

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: PointToOffset_

          PointToOffset_=this%PointToOffset(coord)

          SELECT CASE(this%ID)
          CASE (1)
             SELECT CASE (PointToOffset_)
             CASE (2,4,5,7)
                PointIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (2)
             SELECT CASE (PointToOffset_)
             CASE (1:8)
                PointIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (3)
             SELECT CASE (PointToOffset_)
             CASE (5,11,13,14,16,22)
                PointIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (4)
             SELECT CASE (PointToOffset_)
             CASE (2,4:6,8,10:17,19,21:23,25)
                PointIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (5)
             PointIsInNeighborhood=.TRUE.
             RETURN

          END SELECT

          PointIsInNeighborhood=.FALSE.
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION PointIsInNeighborhood

        LOGICAL FUNCTION OffsetIsInNeighborhood(this,Offset)

          IMPLICIT NONE

          CLASS(Connectivity)    :: this

          INTEGER, INTENT(IN   ) :: Offset

          SELECT CASE(this%ID)
          CASE (1)
             SELECT CASE (Offset)
             CASE (2,4,5,7)
                OffsetIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (2)
             SELECT CASE (Offset)
             CASE (1:8)
                OffsetIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (3)
             SELECT CASE (Offset)
             CASE (5,11,13,14,16,22)
                OffsetIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (4)
             SELECT CASE (Offset)
             CASE (2,4:6,8,10:17,19,21:23,25)
                OffsetIsInNeighborhood=.TRUE.
                RETURN

             END SELECT

          CASE (5)
             OffsetIsInNeighborhood=.TRUE.
             RETURN

          END SELECT

          OffsetIsInNeighborhood=.FALSE.

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION OffsetIsInNeighborhood

        !Convert an offset to a point, in a 3x3x3 cube
        SUBROUTINE OffsetToPoint(this,Offset_,point)

          IMPLICIT NONE

          CLASS(Connectivity)                  :: this

          INTEGER,               INTENT(IN   ) :: Offset_
          INTEGER, DIMENSION(:), INTENT(  OUT) :: Point

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: Offset

          SELECT CASE (this%Dimension)
          CASE (2)
             IF (Offset_.LT.5) THEN
                Offset=Offset_-1
             ELSE
                Offset=Offset_
             ENDIF

             Point(1)=MOD(Offset,3)
             Offset=Offset-Point(1)
             Offset=Offset/3
             Point(1)=Point(1)-1

             Point(2)=MOD(Offset,3)
             Point(2)=Point(2)-1

          CASE (3)
             IF (Offset_.LT.14) THEN
                Offset=Offset_-1
             ELSE
                Offset=Offset_
             ENDIF

             Point(1)=MOD(Offset,3)
             Offset=Offset-Point(1)
             Offset=Offset/3
             Point(1)=Point(1)-1

             Point(2)=MOD(Offset,3)
             Offset=Offset-Point(2)
             Offset=Offset/3
             Point(2)=Point(2)-1

             Point(3)=MOD(Offset,3)
             Point(3)=Point(3)-1

          END SELECT

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE OffsetToPoint

        !Convert a point to an offset, in a 3x3x3 cube
        INTEGER FUNCTION PointToOffset(this,coord)

          IMPLICIT NONE

          CLASS(Connectivity)                  :: this

          INTEGER, DIMENSION(:), INTENT(IN   ) :: coord

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(SIZE(coord)) :: Offset
          INTEGER                         :: info

          CHARACTER(LEN=*), PARAMETER :: caller='PointToOffset'

          Offset=coord+1

          IF (ANY(Offset.GT.2)) THEN
             fail("coords are not in a Unit Cube",ppm_error=ppm_error_fatal,exit_point=no)
          END IF

          SELECT CASE (this%Dimension)
          CASE (2)
             PointToOffset=Offset(1)+3*Offset(2)
             IF (PointToOffset.LT.5) PointToOffset=PointToOffset+1

          CASE (3)
             PointToOffset=Offset(1)+3*(Offset(2)+3*Offset(3))
             IF (PointToOffset.LT.14) PointToOffset=PointToOffset+1

          END SELECT

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION PointToOffset

        FUNCTION FGNeighborhoodConnectivity(this)

          IMPLICIT NONE

          TYPE(Connectivity), POINTER :: this
          TYPE(Connectivity), POINTER :: FGNeighborhoodConnectivity

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller="FGNeighborhoodConnectivity"

          IF (.NOT.ASSOCIATED(this%FGNeighborhood)) THEN
             ALLOCATE(this%FGNeighborhood,STAT=info)
             or_fail_alloc("this%FGNeighborhood",ppm_error=ppm_error_fatal,exit_point=no)

             SELECT CASE (this%Dimension)
             CASE (2)
                ASSOCIATE (FGN => this%FGNeighborhood)
                   CALL FGN%create(2,8,info)
                   or_fail("FGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                END ASSOCIATE
             CASE (3)
                SELECT CASE (this%CellDimension)
                CASE (0)
                   ASSOCIATE (FGN => this%FGNeighborhood)
                      CALL FGN%create(3,26,info)
                      or_fail("FGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                CASE (1,2)
                   ASSOCIATE (FGN => this%FGNeighborhood)
                      CALL FGN%create(3,18,info)
                      or_fail("FGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                END SELECT
             END SELECT !(this%Dimension)
          ENDIF !.NOT.ASSOCIATED(this%FGNeighborhood)

          FGNeighborhoodConnectivity => this%FGNeighborhood

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION FGNeighborhoodConnectivity

        FUNCTION BGNeighborhoodConnectivity(this)

          IMPLICIT NONE

          TYPE(Connectivity), POINTER :: this
          TYPE(Connectivity), POINTER :: BGNeighborhoodConnectivity
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller="BGNeighborhoodConnectivity"

          IF (.NOT.ASSOCIATED(this%BGNeighborhood)) THEN
             ALLOCATE(this%BGNeighborhood,STAT=info)
             or_fail_alloc("this%BGNeighborhood",ppm_error=ppm_error_fatal,exit_point=no)

             SELECT CASE (this%Dimension)
             CASE (2)
                ASSOCIATE (BGN => this%BGNeighborhood)
                   CALL BGN%create(2,8,info)
                   or_fail("BGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                END ASSOCIATE
             CASE (3)
                SELECT CASE (this%CellDimension)
                CASE (0)
                   ASSOCIATE (BGN => this%BGNeighborhood)
                      CALL BGN%create(3,26,info)
                      or_fail("BGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                CASE (1,2)
                   ASSOCIATE (BGN => this%BGNeighborhood)
                      CALL BGN%create(3,18,info)
                      or_fail("BGNeighborhood%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                END SELECT
             END SELECT !(this%Dimension)
          ENDIF !.NOT.ASSOCIATED(this%BGNeighborhood)

          BGNeighborhoodConnectivity => this%BGNeighborhood
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION BGNeighborhoodConnectivity

        FUNCTION BackgroundConnectivity(this)

          IMPLICIT NONE

          TYPE(Connectivity), POINTER :: this
          TYPE(Connectivity), POINTER :: BackgroundConnectivity

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller="BackgroundConnectivity"

          IF (.NOT.ASSOCIATED(this%Background)) THEN
             ALLOCATE(this%Background,STAT=info)
             or_fail_alloc("this%Background",ppm_error=ppm_error_fatal,exit_point=no)

             SELECT CASE (this%Dimension)
             CASE (2)
                SELECT CASE (this%CellDimension)
                CASE (0)
                   ASSOCIATE (BGC => this%Background)
                      CALL BGC%create(2,4,info)
                      or_fail("Background%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                CASE (1)
                   ASSOCIATE (BGC => this%Background)
                      CALL BGC%create(2,8,info)
                      or_fail("Background%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                END SELECT
             CASE (3)
                SELECT CASE (this%CellDimension)
                CASE (0,1)
                   ASSOCIATE (BGC => this%Background)
                      CALL BGC%create(3,6,info)
                      or_fail("Background%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                CASE (2)
                   ASSOCIATE (BGC => this%Background)
                      CALL BGC%create(3,26,info)
                      or_fail("Background%create",ppm_error=ppm_error_fatal,exit_point=no)
                   END ASSOCIATE
                END SELECT
             END SELECT !(this%Dimension)
          ENDIF !.NOT.ASSOCIATED(this%Background)

          BackgroundConnectivity => this%Background
          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION BackgroundConnectivity

        !TOCHECK
        !Using LOGICAL, CONTIGUOUS, DIMENSION(:,:), POINTER :: constructor
        !does not work on intel compiler 13sp1?

        FUNCTION constructor(Connectivity_,NeighborhoodConnectivity_)

          IMPLICIT NONE

          TYPE(Connectivity),      POINTER :: Connectivity_
          TYPE(Connectivity),      POINTER :: NeighborhoodConnectivity_

          LOGICAL, DIMENSION(:,:), POINTER :: constructor

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER               :: nSize,neighbor1,neighbor2
          INTEGER, DIMENSION(2) :: p_2d_1,p_2d_2,sum_2d
          INTEGER, DIMENSION(3) :: p_3d_1,p_3d_2,sum_3d
          INTEGER               :: sumOffset
          INTEGER               :: info

          LOGICAL :: inUnitCube

          CHARACTER(LEN=*), PARAMETER :: caller="UnitCubeNeighbors constructor"

          nSize=Connectivity_%NeighborhoodSize

          ALLOCATE(constructor(nSize,nSize),STAT=info)
          or_fail_alloc("UnitCubeNeighbors allocation Failed!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          constructor=.FALSE.

          IF (Connectivity_%Dimension.NE.NeighborhoodConnectivity_%Dimension) THEN
             fail("check failed: (Connectivity_%Dimension.EQ.NeighborhoodConnectivity_%Dimension) is not true!", &
             & ppm_error=ppm_error_fatal,exit_point=no)
          ENDIF

          SELECT CASE (Connectivity_%Dimension)
          CASE (2)
             DO neighbor1=1,nSize
                IF (NeighborhoodConnectivity_%IsInNeighborhood(neighbor1)) THEN
                   !convert to connectivity OffsetType
                   CALL Connectivity_%OffsetToPoint(neighbor1,p_2d_1)

                   DO neighbor2=1,nSize
                      !convert to Connectivity OffsetType
                      CALL Connectivity_%OffsetToPoint(neighbor2,p_2d_2)

                      sum_2d=p_2d_1+p_2d_2

                      inUnitCube=.NOT.ANY(sum_2d.LT.-1.OR.sum_2d.GT.1)

                      IF (inUnitCube.AND.Connectivity_%AreNeighbors(p_2d_1,sum_2d)) THEN
                         sumOffset=Connectivity_%PointToOffset(sum_2d)
                         constructor(neighbor1,sumOffset)=.TRUE.
                      ENDIF
                   ENDDO !neighbor2

                ENDIF !(NeighborhoodConnectivity_%IsInNeighborhood(neighbor1))
             ENDDO !neighbor1

          CASE (3)
             DO neighbor1=1,nSize
                IF (NeighborhoodConnectivity_%IsInNeighborhood(neighbor1)) THEN
                   !convert to connectivity OffsetType
                   CALL Connectivity_%OffsetToPoint(neighbor1,p_3d_1)

                   DO neighbor2=1,nSize
                      !convert to Connectivity OffsetType
                      CALL Connectivity_%OffsetToPoint(neighbor2,p_3d_2)

                      sum_3d=p_3d_1+p_3d_2

                      inUnitCube=.NOT.ANY(sum_3d.LT.-1.OR.sum_3d.GT.1)

                      IF (inUnitCube.AND.Connectivity_%AreNeighbors(p_3d_1,sum_3d)) THEN
                         sumOffset=Connectivity_%PointToOffset(sum_3d)
                         constructor(neighbor1,sumOffset)=.TRUE.
                      ENDIF
                   ENDDO !neighbor2

                ENDIF !(NeighborhoodConnectivity_%IsInNeighborhood(neighbor1))
             ENDDO !neighbor1

          END SELECT

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION constructor

        SUBROUTINE UnitCubeCCCounter_create(this,Connectivity_,FGNeighborhoodConnectivity_,info)

          IMPLICIT NONE

          CLASS(UnitCubeCCCounter)          :: this

          TYPE(Connectivity), POINTER       :: Connectivity_

          LOGICAL,            INTENT(IN   ) :: FGNeighborhoodConnectivity_

          INTEGER,            INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: i,nSize

          CHARACTER(LEN=*), PARAMETER :: caller="UnitCubeCCCounter_create"

          CALL substart(caller,t0,info)

          this%TConnectivity => Connectivity_

          SELECT CASE (FGNeighborhoodConnectivity_)
          CASE (.TRUE.)
             this%TNeighborhoodConnectivity => FGNeighborhoodConnectivity(Connectivity_)

          CASE DEFAULT
             this%TNeighborhoodConnectivity => BGNeighborhoodConnectivity(Connectivity_)

          END SELECT

          nSize=Connectivity_%NeighborhoodSize

          ALLOCATE(this%ConnectivityTest(nSize),this%m_Image(nSize),STAT=info)
          or_fail_alloc("ConnectivityTest & m_Image")

          DO i=1,nSize
             this%ConnectivityTest(i)=Connectivity_%IsInNeighborhood(i)
          ENDDO

          this%UnitCubeNeighbors => &
          & UnitCubeNeighbors(Connectivity_,this%TNeighborhoodConnectivity)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE UnitCubeCCCounter_create

        SUBROUTINE UnitCubeCCCounter_destroy(this,info)

          IMPLICIT NONE

          CLASS(UnitCubeCCCounter) :: this

          INTEGER,   INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=*), PARAMETER :: caller="UnitCubeCCCounter_destroy"

          CALL substart(caller,t0,info)

          NULLIFY(this%TConnectivity)
          NULLIFY(this%TNeighborhoodConnectivity)

          DEALLOCATE(this%ConnectivityTest,this%m_Image,STAT=info)
          or_fail_dealloc("ConnectivityTest & m_Image")

          DEALLOCATE(this%UnitCubeNeighbors,STAT=info)
          or_fail_dealloc("UnitCubeNeighbors")
          NULLIFY(this%UnitCubeNeighbors)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE UnitCubeCCCounter_destroy

        SUBROUTINE SetImage(this,SubImage,background)

          IMPLICIT NONE

          CLASS(UnitCubeCCCounter)             :: this

          LOGICAL, DIMENSION(:), INTENT(IN   ) :: SubImage

          LOGICAL, OPTIONAL,     INTENT(IN   ) :: background

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: i,SubSize
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller="SetImage"

          SubSize=SIZE(SubImage)

          IF (PRESENT(background)) THEN
             DO i=1,SubSize
                this%m_Image(i)=.NOT.SubImage(i)
             ENDDO
          ELSE
             !!! This check should be done only once, when we are doing the BG
             !!! It is not necessary to do it again
             IF (SIZE(this%m_Image).NE.SubSize) THEN
                fail("(check failed: SIZE(this%m_Image).EQ.SubSize) is not true!", &
                & ppm_error=ppm_error_fatal,exit_point=no)
             ENDIF
             this%m_Image=SubImage
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE SetImage

        INTEGER FUNCTION ConnectedComponentCounter(this)

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list_
          IMPLICIT NONE

          CLASS(UnitCubeCCCounter) :: this

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(ppm_rc_list_) :: seedlst

          INTEGER :: nSize,seed
          INTEGER :: current,neighbor

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: info

          LOGICAL, DIMENSION(:), ALLOCATABLE :: vProcessed_new

          CHARACTER(LEN=*), PARAMETER :: caller="ConnectedComponentCounter"

          nSize=this%TConnectivity%NeighborhoodSize

          !Find first seed
          DO seed=1,nSize
             IF (this%m_Image(seed).AND.this%ConnectivityTest(seed)) EXIT
          ENDDO

          ConnectedComponentCounter=0

          IF (seed.EQ.nSize+1) RETURN

          ALLOCATE(vProcessed_new(nSize),SOURCE=.FALSE.,STAT=info)
          or_fail_alloc("vProcessed_new",ppm_error=ppm_error_fatal,exit_point=no)

          DO WHILE (seed.LE.nSize)
             ConnectedComponentCounter=ConnectedComponentCounter+1
             vProcessed_new(seed)=.TRUE.

             CALL seedlst%add(seed)

             DO WHILE (ASSOCIATED(seedlst%first))
                current = seedlst%first%getValue()
                CALL seedlst%remove()

                DO neighbor=1,nSize
                   IF (.NOT.vProcessed_new(neighbor).AND. &
                   &   this%m_Image(neighbor)       .AND. &
                   &   this%UnitCubeNeighbors(current,neighbor)) THEN
                      CALL seedlst%add(neighbor)
                      vProcessed_new(neighbor)=.TRUE.
                   ENDIF
                ENDDO !neighbor
             ENDDO !(ASSOCIATED(seedlst%first))

             !Look for next seed
             DO seed=1,nSize
                IF (.NOT.vProcessed_new(seed).AND. &
                &   this%m_Image(seed)       .AND. &
                &   this%ConnectivityTest(seed)) EXIT
             ENDDO
          ENDDO !(seed.LE.nSize)

          DEALLOCATE(vProcessed_new,STAT=info)
          or_fail_dealloc("vProcessed_new",ppm_error=ppm_error_fatal,exit_point=no)

          CALL seedlst%destroy()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION ConnectedComponentCounter

        SUBROUTINE ForegroundTopologicalNumberType_create1(this,LabelIn,FGandBGTopoNbPairType_)

          IMPLICIT NONE

          CLASS(ForegroundTopologicalNumberType)     :: this

          INTEGER,                     INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType), INTENT(IN   ) :: FGandBGTopoNbPairType_

          this%label=LabelIn
          this%FGandBGTopoNbPair=FGandBGTopoNbPairType_

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE ForegroundTopologicalNumberType_create1

        SUBROUTINE ForegroundTopologicalNumberType_create2(this,LabelIn,FG,BG)

          IMPLICIT NONE

          CLASS(ForegroundTopologicalNumberType) :: this

          INTEGER,                 INTENT(IN   ) :: LabelIn
          INTEGER,                 INTENT(IN   ) :: FG
          INTEGER,                 INTENT(IN   ) :: BG

          this%label=LabelIn
          this%FGandBGTopoNbPair=FGandBGTopoNbPairType(FG,BG)

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE ForegroundTopologicalNumberType_create2

        SUBROUTINE TopologicalNumberImageFunction_create(this,Connectivity_,info)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          TYPE(Connectivity),     POINTER       :: Connectivity_

          INTEGER,                INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: nSize

          CHARACTER(LEN=*), PARAMETER :: caller="TopologicalNumberImageFunction_create"

          CALL substart(caller,t0,info)

          nSize=Connectivity_%NeighborhoodSize

          ALLOCATE(this%SubImage(nSize),this%DataSubImage(nSize),STAT=info)
          or_fail_alloc("SubImage & DataSubImage")

          this%TFGConnectivity => Connectivity_
          this%TBGConnectivity => BackgroundConnectivity(Connectivity_)

          CALL this%ForegroundUnitCubeCCCounter%create(Connectivity_,.TRUE.,info)
          or_fail("this%ForegroundUnitCubeCCCounter%create")

          CALL this%BackgroundUnitCubeCCCounter%create(this%TBGConnectivity,.FALSE.,info)
          or_fail("this%BackgroundUnitCubeCCCounter%create")

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE TopologicalNumberImageFunction_create

        SUBROUTINE TopologicalNumberImageFunction_destroy(this,info)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER,                INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=*), PARAMETER :: caller="TopologicalNumberImageFunction_destroy"

          CALL substart(caller,t0,info)

          DEALLOCATE(this%SubImage,this%DataSubImage,STAT=info)
          or_fail_dealloc("SubImage & DataSubImage")

          NULLIFY(this%TFGConnectivity)
          NULLIFY(this%TBGConnectivity)

          CALL this%ForegroundUnitCubeCCCounter%destroy(info)
          or_fail("this%ForegroundUnitCubeCCCounter%destroy")

          CALL this%BackgroundUnitCubeCCCounter%destroy(info)
          or_fail("this%BackgroundUnitCubeCCCounter%destroy")

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE TopologicalNumberImageFunction_destroy

        SUBROUTINE readDataSubImagePoint(this,MeshIn,coord)

          USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh_,ppm_t_subpatch_
          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          CLASS(ppm_t_equi_mesh_), POINTER      :: MeshIn

          INTEGER, DIMENSION(:),  INTENT(IN   ) :: coord

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpl_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpl_3d

          INTEGER :: nDim
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller="readDataSubImagePoint"

          nDim=SIZE(coord)

          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             SELECT CASE (nDim)
             CASE (2)
                NULLIFY(wpl_2d)
                !this is the real labels field
                CALL sbpitr%get_field(labels,wpl_2d,info)
                or_fail("Failed to get field data.",ppm_error=ppm_error_fatal,exit_point=no)

                CALL this%readDataSubImage_2d(coord,wpl_2d)

             CASE (3)
                NULLIFY(wpl_3d)
                !this is the real labels field
                CALL sbpitr%get_field(labels,wpl_3d,info)
                or_fail("Failed to get field data.",ppm_error=ppm_error_fatal,exit_point=no)

                CALL this%readDataSubImage_3d(coord,wpl_3d)

             END SELECT
             sbpitr => MeshIn%subpatch%next()
          ENDDO

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE readDataSubImagePoint

        SUBROUTINE readDataSubImage_2d(this,coord,labels_)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)              :: this

          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER :: ll

          ll => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)

          IF (ANY(ll.EQ.FORBIDDEN)) THEN
             IF (ll(1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(1)=0
             ELSE
                this%DataSubImage(1)=ABS(ll(1,1))
             ENDIF
             IF (ll(2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(2)=0
             ELSE
                this%DataSubImage(2)=ABS(ll(2,1))
             ENDIF
             IF (ll(3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(3)=0
             ELSE
                this%DataSubImage(3)=ABS(ll(3,1))
             ENDIF
             IF (ll(1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(4)=0
             ELSE
                this%DataSubImage(4)=ABS(ll(1,2))
             ENDIF
             IF (ll(3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(5)=0
             ELSE
                this%DataSubImage(5)=ABS(ll(3,2))
             ENDIF
             IF (ll(1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(6)=0
             ELSE
                this%DataSubImage(6)=ABS(ll(1,3))
             ENDIF
             IF (ll(2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(7)=0
             ELSE
                this%DataSubImage(7)=ABS(ll(2,3))
             ENDIF
             IF (ll(3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(8)=0
             ELSE
                this%DataSubImage(8)=ABS(ll(3,3))
             ENDIF
          ELSE
             this%DataSubImage(1)=ABS(ll(1,1))
             this%DataSubImage(2)=ABS(ll(2,1))
             this%DataSubImage(3)=ABS(ll(3,1))
             this%DataSubImage(4)=ABS(ll(1,2))
             this%DataSubImage(5)=ABS(ll(3,2))
             this%DataSubImage(6)=ABS(ll(1,3))
             this%DataSubImage(7)=ABS(ll(2,3))
             this%DataSubImage(8)=ABS(ll(3,3))
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE readDataSubImage_2d

        SUBROUTINE readDataSubImage_3d(this,coord,labels_)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)                :: this

          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER :: ll

          ll => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)

          IF (ANY(ll.EQ.FORBIDDEN)) THEN
             IF (ll(1,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(1)=0
             ELSE
                this%DataSubImage(1)=ABS(ll(1,1,1))
             ENDIF
             IF (ll(2,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(2)=0
             ELSE
                this%DataSubImage(2)=ABS(ll(2,1,1))
             ENDIF
             IF (ll(3,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(3)=0
             ELSE
                this%DataSubImage(3)=ABS(ll(3,1,1))
             ENDIF
             IF (ll(1,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(4)=0
             ELSE
                this%DataSubImage(4)=ABS(ll(1,2,1))
             ENDIF
             IF (ll(2,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(5)=0
             ELSE
                this%DataSubImage(5)=ABS(ll(2,2,1))
             ENDIF
             IF (ll(3,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(6)=0
             ELSE
                this%DataSubImage(6)=ABS(ll(3,2,1))
             ENDIF
             IF (ll(1,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(7)=0
             ELSE
                this%DataSubImage(7)=ABS(ll(1,3,1))
             ENDIF
             IF (ll(2,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(8)=0
             ELSE
                this%DataSubImage(8)=ABS(ll(2,3,1))
             ENDIF
             IF (ll(3,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(9)=0
             ELSE
                this%DataSubImage(9)=ABS(ll(3,3,1))
             ENDIF
             IF (ll(1,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(10)=0
             ELSE
                this%DataSubImage(10)=ABS(ll(1,1,2))
             ENDIF
             IF (ll(2,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(11)=0
             ELSE
                this%DataSubImage(11)=ABS(ll(2,1,2))
             ENDIF
             IF (ll(3,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(12)=0
             ELSE
                this%DataSubImage(12)=ABS(ll(3,1,2))
             ENDIF
             IF (ll(1,2,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(13)=0
             ELSE
                this%DataSubImage(13)=ABS(ll(1,2,2))
             ENDIF
             IF (ll(3,2,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(14)=0
             ELSE
                this%DataSubImage(14)=ABS(ll(3,2,2))
             ENDIF
             IF (ll(1,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(15)=0
             ELSE
                this%DataSubImage(15)=ABS(ll(1,3,2))
             ENDIF
             IF (ll(2,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(16)=0
             ELSE
                this%DataSubImage(16)=ABS(ll(2,3,2))
             ENDIF
             IF (ll(3,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(17)=0
             ELSE
                this%DataSubImage(17)=ABS(ll(3,3,2))
             ENDIF
             IF (ll(1,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(18)=0
             ELSE
                this%DataSubImage(18)=ABS(ll(1,1,3))
             ENDIF
             IF (ll(2,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(19)=0
             ELSE
                this%DataSubImage(19)=ABS(ll(2,1,3))
             ENDIF
             IF (ll(3,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(20)=0
             ELSE
                this%DataSubImage(20)=ABS(ll(3,1,3))
             ENDIF
             IF (ll(1,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(21)=0
             ELSE
                this%DataSubImage(21)=ABS(ll(1,2,3))
             ENDIF
             IF (ll(2,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(22)=0
             ELSE
                this%DataSubImage(22)=ABS(ll(2,2,3))
             ENDIF
             IF (ll(3,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(23)=0
             ELSE
                this%DataSubImage(23)=ABS(ll(3,2,3))
             ENDIF
             IF (ll(1,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(24)=0
             ELSE
                this%DataSubImage(24)=ABS(ll(1,3,3))
             ENDIF
             IF (ll(2,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(25)=0
             ELSE
                this%DataSubImage(25)=ABS(ll(2,3,3))
             ENDIF
             IF (ll(3,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(26)=0
             ELSE
                this%DataSubImage(26)=ABS(ll(3,3,3))
             ENDIF
          ELSE
             this%DataSubImage( 1)=ABS(ll(1,1,1))
             this%DataSubImage( 2)=ABS(ll(2,1,1))
             this%DataSubImage( 3)=ABS(ll(3,1,1))
             this%DataSubImage( 4)=ABS(ll(1,2,1))
             this%DataSubImage( 5)=ABS(ll(2,2,1))
             this%DataSubImage( 6)=ABS(ll(3,2,1))
             this%DataSubImage( 7)=ABS(ll(1,3,1))
             this%DataSubImage( 8)=ABS(ll(2,3,1))
             this%DataSubImage( 9)=ABS(ll(3,3,1))
             this%DataSubImage(10)=ABS(ll(1,1,2))
             this%DataSubImage(11)=ABS(ll(2,1,2))
             this%DataSubImage(12)=ABS(ll(3,1,2))
             this%DataSubImage(13)=ABS(ll(1,2,2))
             this%DataSubImage(14)=ABS(ll(3,2,2))
             this%DataSubImage(15)=ABS(ll(1,3,2))
             this%DataSubImage(16)=ABS(ll(2,3,2))
             this%DataSubImage(17)=ABS(ll(3,3,2))
             this%DataSubImage(18)=ABS(ll(1,1,3))
             this%DataSubImage(19)=ABS(ll(2,1,3))
             this%DataSubImage(20)=ABS(ll(3,1,3))
             this%DataSubImage(21)=ABS(ll(1,2,3))
             this%DataSubImage(22)=ABS(ll(2,2,3))
             this%DataSubImage(23)=ABS(ll(3,2,3))
             this%DataSubImage(24)=ABS(ll(1,3,3))
             this%DataSubImage(25)=ABS(ll(2,3,3))
             this%DataSubImage(26)=ABS(ll(3,3,3))
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE readDataSubImage_3d

        SUBROUTINE readDataSubImage__2d(this,tmplabels)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER, DIMENSION(:,:), POINTER      :: tmplabels

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          IF (ANY(tmplabels.EQ.FORBIDDEN)) THEN
             IF (tmplabels(1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(1)=0
             ELSE
                this%DataSubImage(1)=ABS(tmplabels(1,1))
             ENDIF
             IF (tmplabels(2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(2)=0
             ELSE
                this%DataSubImage(2)=ABS(tmplabels(2,1))
             ENDIF
             IF (tmplabels(3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(3)=0
             ELSE
                this%DataSubImage(3)=ABS(tmplabels(3,1))
             ENDIF
             IF (tmplabels(1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(4)=0
             ELSE
                this%DataSubImage(4)=ABS(tmplabels(1,2))
             ENDIF
             IF (tmplabels(3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(5)=0
             ELSE
                this%DataSubImage(5)=ABS(tmplabels(3,2))
             ENDIF
             IF (tmplabels(1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(6)=0
             ELSE
                this%DataSubImage(6)=ABS(tmplabels(1,3))
             ENDIF
             IF (tmplabels(2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(7)=0
             ELSE
                this%DataSubImage(7)=ABS(tmplabels(2,3))
             ENDIF
             IF (tmplabels(3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(8)=0
             ELSE
                this%DataSubImage(8)=ABS(tmplabels(3,3))
             ENDIF
          ELSE
             this%DataSubImage(1)=ABS(tmplabels(1,1))
             this%DataSubImage(2)=ABS(tmplabels(2,1))
             this%DataSubImage(3)=ABS(tmplabels(3,1))
             this%DataSubImage(4)=ABS(tmplabels(1,2))
             this%DataSubImage(5)=ABS(tmplabels(3,2))
             this%DataSubImage(6)=ABS(tmplabels(1,3))
             this%DataSubImage(7)=ABS(tmplabels(2,3))
             this%DataSubImage(8)=ABS(tmplabels(3,3))
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE readDataSubImage__2d

        SUBROUTINE readDataSubImage__3d(this,tmplabels)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER, DIMENSION(:,:,:), POINTER    :: tmplabels
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          IF (ANY(tmplabels.EQ.FORBIDDEN)) THEN
             IF (tmplabels(1,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(1)=0
             ELSE
                this%DataSubImage(1)=ABS(tmplabels(1,1,1))
             ENDIF
             IF (tmplabels(2,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(2)=0
             ELSE
                this%DataSubImage(2)=ABS(tmplabels(2,1,1))
             ENDIF
             IF (tmplabels(3,1,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(3)=0
             ELSE
                this%DataSubImage(3)=ABS(tmplabels(3,1,1))
             ENDIF
             IF (tmplabels(1,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(4)=0
             ELSE
                this%DataSubImage(4)=ABS(tmplabels(1,2,1))
             ENDIF
             IF (tmplabels(2,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(5)=0
             ELSE
                this%DataSubImage(5)=ABS(tmplabels(2,2,1))
             ENDIF
             IF (tmplabels(3,2,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(6)=0
             ELSE
                this%DataSubImage(6)=ABS(tmplabels(3,2,1))
             ENDIF
             IF (tmplabels(1,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(7)=0
             ELSE
                this%DataSubImage(7)=ABS(tmplabels(1,3,1))
             ENDIF
             IF (tmplabels(2,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(8)=0
             ELSE
                this%DataSubImage(8)=ABS(tmplabels(2,3,1))
             ENDIF
             IF (tmplabels(3,3,1).EQ.FORBIDDEN) THEN
                this%DataSubImage(9)=0
             ELSE
                this%DataSubImage(9)=ABS(tmplabels(3,3,1))
             ENDIF
             IF (tmplabels(1,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(10)=0
             ELSE
                this%DataSubImage(10)=ABS(tmplabels(1,1,2))
             ENDIF
             IF (tmplabels(2,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(11)=0
             ELSE
                this%DataSubImage(11)=ABS(tmplabels(2,1,2))
             ENDIF
             IF (tmplabels(3,1,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(12)=0
             ELSE
                this%DataSubImage(12)=ABS(tmplabels(3,1,2))
             ENDIF
             IF (tmplabels(1,2,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(13)=0
             ELSE
                this%DataSubImage(13)=ABS(tmplabels(1,2,2))
             ENDIF
             IF (tmplabels(3,2,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(14)=0
             ELSE
                this%DataSubImage(14)=ABS(tmplabels(3,2,2))
             ENDIF
             IF (tmplabels(1,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(15)=0
             ELSE
                this%DataSubImage(15)=ABS(tmplabels(1,3,2))
             ENDIF
             IF (tmplabels(2,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(16)=0
             ELSE
                this%DataSubImage(16)=ABS(tmplabels(2,3,2))
             ENDIF
             IF (tmplabels(3,3,2).EQ.FORBIDDEN) THEN
                this%DataSubImage(17)=0
             ELSE
                this%DataSubImage(17)=ABS(tmplabels(3,3,2))
             ENDIF
             IF (tmplabels(1,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(18)=0
             ELSE
                this%DataSubImage(18)=ABS(tmplabels(1,1,3))
             ENDIF
             IF (tmplabels(2,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(19)=0
             ELSE
                this%DataSubImage(19)=ABS(tmplabels(2,1,3))
             ENDIF
             IF (tmplabels(3,1,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(20)=0
             ELSE
                this%DataSubImage(20)=ABS(tmplabels(3,1,3))
             ENDIF
             IF (tmplabels(1,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(21)=0
             ELSE
                this%DataSubImage(21)=ABS(tmplabels(1,2,3))
             ENDIF
             IF (tmplabels(2,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(22)=0
             ELSE
                this%DataSubImage(22)=ABS(tmplabels(2,2,3))
             ENDIF
             IF (tmplabels(3,2,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(23)=0
             ELSE
                this%DataSubImage(23)=ABS(tmplabels(3,2,3))
             ENDIF
             IF (tmplabels(1,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(24)=0
             ELSE
                this%DataSubImage(24)=ABS(tmplabels(1,3,3))
             ENDIF
             IF (tmplabels(2,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(25)=0
             ELSE
                this%DataSubImage(25)=ABS(tmplabels(2,3,3))
             ENDIF
             IF (tmplabels(3,3,3).EQ.FORBIDDEN) THEN
                this%DataSubImage(26)=0
             ELSE
                this%DataSubImage(26)=ABS(tmplabels(3,3,3))
             ENDIF
          ELSE
             this%DataSubImage( 1)=ABS(tmplabels(1,1,1))
             this%DataSubImage( 2)=ABS(tmplabels(2,1,1))
             this%DataSubImage( 3)=ABS(tmplabels(3,1,1))
             this%DataSubImage( 4)=ABS(tmplabels(1,2,1))
             this%DataSubImage( 5)=ABS(tmplabels(2,2,1))
             this%DataSubImage( 6)=ABS(tmplabels(3,2,1))
             this%DataSubImage( 7)=ABS(tmplabels(1,3,1))
             this%DataSubImage( 8)=ABS(tmplabels(2,3,1))
             this%DataSubImage( 9)=ABS(tmplabels(3,3,1))
             this%DataSubImage(10)=ABS(tmplabels(1,1,2))
             this%DataSubImage(11)=ABS(tmplabels(2,1,2))
             this%DataSubImage(12)=ABS(tmplabels(3,1,2))
             this%DataSubImage(13)=ABS(tmplabels(1,2,2))
             this%DataSubImage(14)=ABS(tmplabels(3,2,2))
             this%DataSubImage(15)=ABS(tmplabels(1,3,2))
             this%DataSubImage(16)=ABS(tmplabels(2,3,2))
             this%DataSubImage(17)=ABS(tmplabels(3,3,2))
             this%DataSubImage(18)=ABS(tmplabels(1,1,3))
             this%DataSubImage(19)=ABS(tmplabels(2,1,3))
             this%DataSubImage(20)=ABS(tmplabels(3,1,3))
             this%DataSubImage(21)=ABS(tmplabels(1,2,3))
             this%DataSubImage(22)=ABS(tmplabels(2,2,3))
             this%DataSubImage(23)=ABS(tmplabels(3,2,3))
             this%DataSubImage(24)=ABS(tmplabels(1,3,3))
             this%DataSubImage(25)=ABS(tmplabels(2,3,3))
             this%DataSubImage(26)=ABS(tmplabels(3,3,3))
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END SUBROUTINE readDataSubImage__3d

        FUNCTION EvaluateAtIndex1_2d(this,coord,labels_) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)              :: this

          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_

          TYPE(FGandBGTopoNbPairType)                        :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(coord,labels_)

!           this%SubImage=MERGE(.FALSE.,.TRUE.,this%DataSubImage.EQ.0)
          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex1_2d

        FUNCTION EvaluateAtIndex1_3d(this,coord,labels_) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)                :: this

          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_

          TYPE(FGandBGTopoNbPairType)                          :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(coord,labels_)

          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex1_3d

        FUNCTION EvaluateAtIndex1__2d(this,tmplabels) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER, DIMENSION(:,:), POINTER      :: tmplabels

          TYPE(FGandBGTopoNbPairType)           :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(tmplabels)

          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex1__2d

        FUNCTION EvaluateAtIndex1__3d(this,tmplabels) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER, DIMENSION(:,:,:), POINTER    :: tmplabels

          TYPE(FGandBGTopoNbPairType)           :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(tmplabels)

          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex1__3d

        FUNCTION EvaluateAtIndex2(this,MeshIn,coord) RESULT(FGBGTNP)

          USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh_
          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)  :: this

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn

          INTEGER, DIMENSION(:),   INTENT(IN   ) :: coord

          TYPE(FGandBGTopoNbPairType)            :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(MeshIn,coord)

          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex2

        FUNCTION EvaluateAtIndex3(this) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          TYPE(FGandBGTopoNbPairType)           :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          this%SubImage=this%DataSubImage.NE.0

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          IF (this%ComputeForegroundTN) THEN
             FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%FGNumber = 0
          ENDIF

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          IF (this%ComputeBackgroundTN) THEN
             FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()
          ELSE
             FGBGTNP%BGNumber = 0
          ENDIF

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAtIndex3

        FUNCTION EvaluateFGTNOfLabelAtIndex1_2d(this,coord,labels_,LabelIn) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)              :: this

          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,                             INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)                        :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(coord,labels_)

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex1_2d

        FUNCTION EvaluateFGTNOfLabelAtIndex1_3d(this,coord,labels_,LabelIn) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)                :: this

          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,                               INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)                          :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(coord,labels_)

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex1_3d

        FUNCTION EvaluateFGTNOfLabelAtIndex1__2d(this,tmplabels,LabelIn) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)  :: this

          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)            :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(tmplabels)

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex1__2d

        FUNCTION EvaluateFGTNOfLabelAtIndex1__3d(this,tmplabels,LabelIn) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)    :: this

          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)              :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(tmplabels)

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex1__3d

        FUNCTION EvaluateFGTNOfLabelAtIndex2(this,MeshIn,coord,LabelIn) RESULT(FGBGTNP)

          USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh_
          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)  :: this

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn

          INTEGER, DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,                 INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)            :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          CALL this%readDataSubImage(MeshIn,coord)

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber = this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber = this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex2

        FUNCTION EvaluateFGTNOfLabelAtIndex3(this,LabelIn) RESULT(FGBGTNP)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction) :: this

          INTEGER,                INTENT(IN   ) :: LabelIn

          TYPE(FGandBGTopoNbPairType)           :: FGBGTNP

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          this%SubImage=this%DataSubImage.EQ.LabelIn

          CALL this%ForegroundUnitCubeCCCounter%SetImage(this%SubImage)

          FGBGTNP%FGNumber=this%ForegroundUnitCubeCCCounter%ConnectedComponentCounter()

          !Topological number in the background
          CALL this%BackgroundUnitCubeCCCounter%SetImage(this%SubImage,background=.TRUE.)

          FGBGTNP%BGNumber=this%BackgroundUnitCubeCCCounter%ConnectedComponentCounter()

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateFGTNOfLabelAtIndex3

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1_2d(this,coord,labels_) RESULT(vTNvector)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)              :: this

          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          CALL this%readDataSubImage(coord,labels_)

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1_2d

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1_3d(this,coord,labels_) RESULT(vTNvector)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)                :: this

          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          CALL this%readDataSubImage(coord,labels_)

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1_3d

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1__2d(this,tmplabels) RESULT(vTNvector)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)        :: this

          INTEGER,             DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          CALL this%readDataSubImage(tmplabels)

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1__2d

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1__3d(this,tmplabels) RESULT(vTNvector)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)          :: this

          INTEGER,             DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          CALL this%readDataSubImage(tmplabels)

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex1__3d

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex2(this,MeshIn,coord) RESULT(vTNvector)

          USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh_
          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)                  :: this

          CLASS(ppm_t_equi_mesh_), POINTER                       :: MeshIn

          INTEGER,                 DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS,     DIMENSION(:,:), POINTER       :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          CALL this%readDataSubImage(MeshIn,coord)

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex2

        FUNCTION EvaluateAdjacentRegionsFGTNAtIndex3(this) RESULT(vTNvector)

          IMPLICIT NONE

          CLASS(TopologicalNumberImageFunction)        :: this

          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER :: vTNvector

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP


          INTEGER, DIMENSION(26) :: vAdjacentLabels
          INTEGER                :: nSize,i,j,vLinearIndex,vLabeln
          INTEGER                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="EvaluateAdjacentRegionsFGTNAtIndex"

          nSize=this%TFGConnectivity%NumberOfNeighbors

          j=0
          DO i=1,nSize
             vLinearIndex=this%TFGConnectivity%NeighborsOffsets(i)
             IF (this%DataSubImage(vLinearIndex).NE.0) THEN
                j=j+1
                vLabeln=this%DataSubImage(vLinearIndex)
                vAdjacentLabels(j)=vLabeln
                IF (ANY(vLabeln.EQ.vAdjacentLabels(1:j-1))) j=j-1
             ENDIF
          ENDDO

          ALLOCATE(vTNvector(3,j),STAT=info)
          or_fail_alloc("vTNvector",ppm_error=ppm_error_fatal,exit_point=no)

          DO i=1,j
             vLabeln = vAdjacentLabels(i)

             FGBGTNP=this%EvaluateFGTNOfLabelAtIndex(vLabeln)

             vTNvector(1,i)=vLabeln
             vTNvector(2,i)=FGBGTNP%FGNumber
             vTNvector(3,i)=FGBGTNP%BGNumber
          ENDDO !

          !------------------------------------------------------------------
          !  Return
          !------------------------------------------------------------------
          RETURN
        END FUNCTION EvaluateAdjacentRegionsFGTNAtIndex3

        LOGICAL FUNCTION IsSingleFGPoint_2d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,                             INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsSingleFGPoint_2d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2))).EQ.LabelIn) THEN
                IsSingleFGPoint_2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsSingleFGPoint_2d

        LOGICAL FUNCTION IsSingleFGPoint__2d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsSingleFGPoint__2d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2))).EQ.LabelIn) THEN
                IsSingleFGPoint__2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsSingleFGPoint__2d

        LOGICAL FUNCTION IsSingleFGPoint_3d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,                               INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsSingleFGPoint_3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2),ll(3))).EQ.LabelIn) THEN
                IsSingleFGPoint_3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsSingleFGPoint_3d

        LOGICAL FUNCTION IsSingleFGPoint__3d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsSingleFGPoint__3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2),ll(3))).EQ.LabelIn) THEN
                IsSingleFGPoint__3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsSingleFGPoint__3d

        LOGICAL FUNCTION IsEnclosedByLabel_FGConnectivity_2d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,                             INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_FGConnectivity_2d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2))).NE.LabelIn) THEN
                IsEnclosedByLabel_FGConnectivity_2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_FGConnectivity_2d

        LOGICAL FUNCTION IsEnclosedByLabel_FGConnectivity__2d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_FGConnectivity__2d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2))).NE.LabelIn) THEN
                IsEnclosedByLabel_FGConnectivity__2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_FGConnectivity__2d

        LOGICAL FUNCTION IsEnclosedByLabel_FGConnectivity_3d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,                               INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_FGConnectivity_3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2),ll(3))).NE.LabelIn) THEN
                IsEnclosedByLabel_FGConnectivity_3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_FGConnectivity_3d

        LOGICAL FUNCTION IsEnclosedByLabel_FGConnectivity__3d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_FGConnectivity__3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2),ll(3))).NE.LabelIn) THEN
                IsEnclosedByLabel_FGConnectivity__3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_FGConnectivity__3d

        LOGICAL FUNCTION IsEnclosedByLabel_BGConnectivity_2d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,                             INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_BGConnectivity_2d=.TRUE.
          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=coord+BG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2))).NE.LabelIn) THEN
                IsEnclosedByLabel_BGConnectivity_2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_BGConnectivity_2d

        LOGICAL FUNCTION IsEnclosedByLabel_BGConnectivity__2d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_BGConnectivity__2d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2))).NE.LabelIn) THEN
                IsEnclosedByLabel_BGConnectivity__2d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_BGConnectivity__2d

        LOGICAL FUNCTION IsEnclosedByLabel_BGConnectivity_3d(coord,labels_,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,                               INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_BGConnectivity_3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
             IF (ABS(labels_(ll(1),ll(2),ll(3))).NE.LabelIn) THEN
                IsEnclosedByLabel_BGConnectivity_3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_BGConnectivity_3d

        LOGICAL FUNCTION IsEnclosedByLabel_BGConnectivity__3d(tmplabels,LabelIn)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: LabelIn
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: i

          IsEnclosedByLabel_BGConnectivity__3d=.TRUE.
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
             IF (ABS(tmplabels(ll(1),ll(2),ll(3))).NE.LabelIn) THEN
                IsEnclosedByLabel_BGConnectivity__3d=.FALSE.
                RETURN
             ENDIf
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION IsEnclosedByLabel_BGConnectivity__3d

