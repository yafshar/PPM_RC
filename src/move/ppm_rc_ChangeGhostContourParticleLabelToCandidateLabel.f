      !-------------------------------------------------------------------------
      !  Subroutine   :   ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel
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
      !  Subroutine   :   ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Loops through all the particle Candidates, and remove
      !                 particles which are illegal
      !
      !                The method performs all steps needed to switch a pixel to
      !                a different label. No tests about topology or validity of
      !                the arguments are performed.
      !                Concretely,
      !                - the label image is updated.
      !                - the main contour container is maintained.
      !                - Statistics are updated, including energy specific tasks
      !                  (maintaining the deconvolution-model image.
      !                aNewCCPointsContainer contains contour points that have
      !                been added to the
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel) &
      &          (ipatch,labels_,pind_,coord,Nm,vFromLabel,vToLabel,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,                                INTENT(IN   ) :: ipatch
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
        INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: pind_
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: pind_
#endif
        INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
        INTEGER,             DIMENSION(:),      POINTER       :: Nm
        INTEGER,                                INTENT(IN   ) :: vFromLabel
        INTEGER,                                INTENT(IN   ) :: vToLabel
        INTEGER,                                INTENT(  OUT) :: info
        !-----------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        TYPE(ppm_rc_list), POINTER :: seed

        INTEGER, DIMENSION(__DIME) :: ll,kk
        INTEGER               :: i,j,qpart,val
        INTEGER               :: vLabel
        INTEGER               :: nsize

        LOGICAL :: IsEnclosedByLabelFGConnectivity

        start_subroutine("ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel")

        IF (vFromLabel.EQ.vToLabel) GOTO 9999

        IF (vFromLabel.NE.0) THEN
           !It would be Greater than zero
           ! The loser region (if it is not the BG region) has to add the
           ! neighbors of the lost point to the contour list.
           !
           !adding all the necessary neighbors to the contour container
           !when removing a point from the region and adding it to the background
           !region.
           !Basically all the FG-Neighborhood neighbors with the same label are added
           !to the container.
           DO i=1,FG_ConnectivityType%NumberOfNeighbors
              ll=coord+FG_ConnectivityType%NeighborsPoints(:,i)
              !TOCHECK
              !TODO
              IF (ANY(ll.LT.1.OR.ll.GT.Nm)) CYCLE
#if   __DIME == __2D
              vLabel=labels_(ll(1),ll(2))
#elif __DIME == __3D
              vLabel=labels_(ll(1),ll(2),ll(3))
#endif

              IF (vLabel.EQ.vFromLabel) THEN
                 !New point is an inner point with the same label as the contour
                 !particle
                 NpartNewtmp=NpartNewtmp+1
                 qpart=Mpart+NpartNew+NpartNewtmp
#if   __DIME == __2D
                 labels_(ll(1),ll(2))=-vFromLabel
                 pind_(ll(1),ll(2))=qpart
#elif __DIME == __3D
                 labels_(ll(1),ll(2),ll(3))=-vFromLabel
                 pind_(ll(1),ll(2),ll(3))=qpart
#endif

                 ALLOCATE(seed,STAT=info)
                 or_fail_alloc("seed")
                 CALL seed%add(ll)
                 CALL InnerContourContainer(ipatch)%push(seed,info)
                 or_fail("could not add new seed to the collection")

                 IF (NpartNew+NpartNewtmp.GT.htableg%nrow) THEN
                    CALL htableg%grow(info)
                    or_fail("htableg%grow")
                 ENDIF

                 CALL htableg%insert(qpart,InnerContourContainer(ipatch)%max_id,info)
                 or_fail("hash insert")
              ENDIF
           ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
        ENDIF !(vFromLabel.GT.0)

        !
        ! Erase the point from the surface container in case it now belongs to
        ! the background. Else, add the point to the container (or replace it
        ! in case it has been there already).
        !
        IF (vToLabel.EQ.0) THEN
#if   __DIME == __2D
           pind_(coord(1),coord(2))=0
#elif __DIME == __3D
           pind_(coord(1),coord(2),coord(3))=0
#endif
        ELSE
           nsize=SIZE(list_del_parts)

           ! Remove 'enclosed' contour points from the container.
           ! For the BG this makes no sense.
           !
           ! Maintain the inner contour container:
           ! - Remove all the indices in the BG-connected neighborhood, that are
           !   interior points, from the contour container.
           !   Interior here means that none neighbors in the FG-Neighborhood
           !   has a different label.
           !
           DO i=1,BG_ConnectivityType%NumberOfNeighbors
              ll=coord+BG_ConnectivityType%NeighborsPoints(:,i)
              IF (ANY(ll.LT.1.OR.ll.GT.Nm)) CYCLE
#if   __DIME == __2D
              vLabel=labels_(ll(1),ll(2))
#elif __DIME == __3D
              vLabel=labels_(ll(1),ll(2),ll(3))
#endif
              IF (vLabel.EQ.-vToLabel) THEN
                 IsEnclosedByLabelFGConnectivity=.TRUE.
                 DO j=1,FG_ConnectivityType%NumberOfNeighbors
                    kk=ll+FG_ConnectivityType%NeighborsPoints(:,j)
#if   __DIME == __2D
                    IF (ABS(labels_(kk(1),kk(2))).NE.vToLabel) THEN
#elif __DIME == __3D
                    IF (ABS(labels_(kk(1),kk(2),kk(3))).NE.vToLabel) THEN
#endif
                       IsEnclosedByLabelFGConnectivity=.FALSE.
                       EXIT
                    ENDIF
                 ENDDO !j=1,FG_ConnectivityType%NumberOfNeighbors
                 IF (IsEnclosedByLabelFGConnectivity) THEN
#if   __DIME == __2D
                    labels_(ll(1),ll(2))=vToLabel
                    qpart=pind_(ll(1),ll(2))
                    pind_(ll(1),ll(2))=0
#elif __DIME == __3D
                    labels_(ll(1),ll(2),ll(3))=vToLabel
                    qpart=pind_(ll(1),ll(2),ll(3))
                    pind_(ll(1),ll(2),ll(3))=0
#endif
                    IF (qpart.GT.0.AND.qpart.LE.Npart) THEN
                       del_parts=del_parts+1
                       IF (del_parts.GT.nsize) THEN
                          ALLOCATE(bufi(nsize*2),STAT=info)
                          or_fail_alloc("Tmp array allocation failed!")

                          FORALL (j=1:nsize) bufi(j)=list_del_parts(j)

                          nsize=nsize*2

                          CALL MOVE_ALLOC(bufi,list_del_parts)
                       ENDIF
                       list_del_parts(del_parts)=qpart
                    ELSE IF (qpart.GT.Mpart) THEN
                       val=htableg%search(qpart)
                       IF (val.NE.htable_null) THEN
                          del_parts=del_parts+1
                          IF (del_parts.GT.nsize) THEN
                             ALLOCATE(bufi(nsize*2),STAT=info)
                             or_fail_alloc("Tmp array allocation failed!")

                             FORALL (j=1:nsize) bufi(j)=list_del_parts(j)

                             nsize=nsize*2

                             CALL MOVE_ALLOC(bufi,list_del_parts)
                          ENDIF
                          list_del_parts(del_parts)=val
                       ENDIF !val.NE.htable_null
                    ENDIF !qpart.GT.0.AND.qpart.LE.Npart
                 ENDIF !(IsEnclosedByLabelFGConnectivity)
              ENDIF !(vLabel.EQ.-vToLabel)
           ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
        ENDIF !(vToLabel.EQ.0)

        end_subroutine()

      END SUBROUTINE DTYPE(ppm_rc_ChangeGhostContourParticleLabelToCandidateLabel)
