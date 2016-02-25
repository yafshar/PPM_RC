      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_floodfill
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
      !  Author           - y.afshar           April   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_floodfill
      !-------------------------------------------------------------------------
      !
      !  Purpose      : The purpose of floodFill is to color an entire area
      !                 of connected pixels with the same color
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       :
      !                 info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_floodFillScanline)(labelled, &
      &          labels_,Nm_,coord,oldlabel,newlabel,Sm,info)
        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:,:),   ALLOCATABLE   :: labelled
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: labels_
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
#endif
        INTEGER, DIMENSION(:),     POINTER       :: Nm_
        !!! number of mesh points
        INTEGER, DIMENSION(:),     INTENT(IN   ) :: coord
        INTEGER,                   INTENT(IN   ) :: oldlabel
        INTEGER,                   INTENT(IN   ) :: newlabel
        INTEGER,                   INTENT(IN   ) :: Sm
        !!! mesh starting point number whether it will start from 0(ghost) or 1
        INTEGER,                   INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_rc_list) :: seedlst
        TYPE(ppm_rc_list) :: seedlsti

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(:,:), ALLOCATABLE :: tlabelled
        INTEGER, DIMENSION(:),   POINTER     :: seedn
        INTEGER, DIMENSION(__DIME)                :: ld
#if   __DIME == __3D
        INTEGER, DIMENSION(__DIME)                :: ldi
#endif
        INTEGER, DIMENSION(__DIME)                :: Nm,Nmm
        INTEGER                              :: xx,xs,xe
        INTEGER                              :: nsize,nsizet,nb,Smp
        INTEGER                              :: ii,jj

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_floodFillScanline"

        LOGICAL :: south,north
#if   __DIME == __3D
        LOGICAL :: top,bottom
#endif
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        nb=0
        !!! number of pixels(or voxels) in a cluster

        IF (newlabel.EQ.oldlabel) THEN
           IF (ALLOCATED(labelled)) THEN
              nsize=SIZE(labelled,DIM=2)
           ELSE
              nsize=8
              ALLOCATE(labelled(__DIME,nsize),STAT=info)
              or_fail_alloc("labelled")
           ENDIF
           GOTO 8888
        ENDIF

        Nm=Nm_(1:__DIME)+MERGE(1,0,Sm.EQ.0)
        Nmm=Nm-1

        Smp=Sm+1

#if   __DIME == __2D
        IF (ABS(labels_(coord(1),coord(2))).EQ.oldlabel) THEN
#elif __DIME == __3D
        IF (ABS(labels_(coord(1),coord(2),coord(3))).EQ.oldlabel) THEN
#endif
           IF (ALLOCATED(labelled)) THEN
              nsize=SIZE(labelled,DIM=2)
           ELSE
              nsize=64
              ALLOCATE(labelled(__DIME,nsize),STAT=info)
              or_fail_alloc("labelled")
           ENDIF
#if   __DIME == __2D
        ELSE IF (ABS(labels_(coord(1),coord(2))).EQ.newlabel) THEN
#elif __DIME == __3D
        ELSE IF (ABS(labels_(coord(1),coord(2),coord(3))).EQ.newlabel) THEN
#endif
           IF (ALLOCATED(labelled)) THEN
              DEALLOCATE(labelled,STAT=info)
              or_fail_dealloc("labelled")
           ENDIF
           nsize=16
           ALLOCATE(labelled(__DIME,nsize),STAT=info)
           or_fail_alloc("labelled")
#if   __DIME == __2D
           labels_(coord(1),coord(2))=SIGN(oldlabel,labels_(coord(1),coord(2)))
#elif __DIME == __3D
           labels_(coord(1),coord(2),coord(3))=SIGN(oldlabel,labels_(coord(1),coord(2),coord(3)))
#endif
        ENDIF

        !TODO
        !it only has been implemented for 6 and 4 connectivity forground
        !connectivity(3DIM and 2DIM respectively.)
        CALL seedlst%add(coord)
        SELECT CASE (FG_ConnectivityType%ID)
        CASE (1,3)
           seedlst_loop: DO WHILE (ASSOCIATED(seedlst%first))
              seedn => seedlst%first%getValue()
              ld=seedn
              CALL seedlst%remove()

              IF (ANY(ld.LT.Sm.OR.ld.GT.Nm)) THEN
                 CYCLE seedlst_loop
              ENDIF

#if   __DIME == __2D
              DO xs=ld(1),Sm,-1
                 IF (ABS(labels_(xs,ld(2))).NE.oldlabel) EXIT
              ENDDO
              xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.ABS(labels_(Sm,ld(2))).EQ.oldlabel)

              DO xe=ld(1)+1,Nm(1)
                 IF (ABS(labels_(xe,ld(2))).NE.oldlabel) EXIT
              ENDDO
              xe=xe-1

              IF (nsize.LT.nb+xe-xs+1) THEN
                 nsizet=nsize
                 DO WHILE (nsizet.LT.nb+xe-xs+1)
                    nsizet=nsizet*2
                 ENDDO
                 ALLOCATE(tlabelled(2,nsizet),STAT=info)
                 or_fail_alloc("tlabelled")

                 FORALL (ii=1:2,jj=1:nsize) tlabelled(ii,jj)=labelled(ii,jj)

                 nsize=nsizet

                 CALL MOVE_ALLOC(tlabelled,labelled)
              ENDIF
              south=.FALSE.
              north=.FALSE.
              IF      (ld(2).EQ.Sm) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))
                    nb=nb+1
                    labelled(1,nb)=xx
                    labelled(2,nb)=ld(2)

                    IF (.NOT.north.AND.ABS(labels_(xx,Smp)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,Smp)
                       north=.TRUE.
                    ELSE IF (north.AND.ABS(labels_(xx,Smp)).NE.oldlabel) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE IF (ld(2).EQ.Nm(2)) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))
                    nb=nb+1
                    labelled(1,nb)=xx
                    labelled(2,nb)=ld(2)

                    IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2))).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,Nmm(2))
                       south=.TRUE.
                    ELSE IF (south.AND.ABS(labels_(xx,Nmm(2))).NE.oldlabel) THEN
                       south=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))
                    nb=nb+1
                    labelled(1,nb)=xx
                    labelled(2,nb)=ld(2)

                    IF (.NOT.south.AND.ABS(labels_(xx,ld(2)-1)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,ld(2)-1)
                       south=.TRUE.
                    ELSE IF (south.AND.ABS(labels_(xx,ld(2)-1)).NE.oldlabel) THEN
                       south=.FALSE.
                    ENDIF
                    IF (.NOT.north.AND.ABS(labels_(xx,ld(2)+1)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,ld(2)+1)
                       north=.TRUE.
                    ELSE IF (north.AND.ABS(labels_(xx,ld(2)+1)).NE.oldlabel) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ENDIF
#elif __DIME == __3D
              CALL seedlsti%add(ld)
              seedlsti_loop: DO WHILE (ASSOCIATED(seedlsti%first))
                 seedn => seedlsti%first%getValue()
                 ldi=seedn
                 CALL seedlsti%remove()

                 IF (ANY(ldi.LT.Sm.OR.ldi.GT.Nm)) THEN
                    CYCLE seedlsti_loop
                 ENDIF

                 DO xs=ldi(1),Sm,-1
                    IF (ABS(labels_(xs,ldi(2),ldi(3))).NE.oldlabel) EXIT
                 ENDDO
                 xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.ABS(labels_(Sm,ldi(2),ldi(3))).EQ.oldlabel)

                 DO xe=ldi(1)+1,Nm(1)
                    IF (ABS(labels_(xe,ldi(2),ldi(3))).NE.oldlabel) EXIT
                 ENDDO
                 xe=xe-1

                 IF (nsize.LT.nb+xe-xs+1) THEN
                    nsizet=nsize
                    DO WHILE (nsizet.LT.nb+xe-xs+1)
                       nsizet=nsizet*2
                    ENDDO
                    ALLOCATE(tlabelled(3,nsizet),STAT=info)
                    or_fail_alloc("tlabelled")

                    FORALL (ii=1:3,jj=1:nsize) tlabelled(ii,jj)=labelled(ii,jj)

                    nsize=nsizet

                    CALL MOVE_ALLOC(tlabelled,labelled)
                 ENDIF

                 south =.FALSE.
                 north =.FALSE.
                 bottom=.FALSE.
                 top   =.FALSE.
                 IF      (ldi(3).EQ.Sm) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,Sm)).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,Sm)
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,Sm)).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,Sm,Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,Sm,Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,Sm,Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !(ldi(2).EQ.Sm)
                 ELSE IF (ldi(3).EQ.Nm(3)) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ELSE
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))
                          nb=nb+1
                          labelled(1,nb)=xx
                          labelled(2:3,nb)=ldi(2:3)

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ENDIF !ldi(3).EQ.Sm
              ENDDO seedlsti_loop !ASSOCIATED(seedlsti%first)
#endif
           ENDDO seedlst_loop !ASSOCIATED(seedlst%first)

        CASE DEFAULT
           fail("It has not implemented for this Foreground connectivity",ppm_error=ppm_error_fatal)

        END SELECT

        CALL seedlst%destroy()
        CALL seedlsti%destroy()

      8888 CONTINUE

        nb=nb+1
        IF (nsize.GE.nb) THEN
           FORALL (ii=1:__DIME,jj=nb:nsize) labelled(ii,jj)=-1
           !coordinate in labelled are between zero and Nm+1
           !so the -1 value can seperate the non labelled ones
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_floodFillScanline)

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_floodfill
      !-------------------------------------------------------------------------
      !
      !  Purpose      : The purpose of floodFill is to color an entire area
      !                 of connected pixels with the same color
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       :
      !                 info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Author           - y.afshar           April   2014
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_floodFillScanline_)(labels_, &
      &          Nm_,coord,oldlabel,newlabel,Sm,info)
        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: labels_
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
#endif
        INTEGER, DIMENSION(:),     POINTER       :: Nm_
        !!! number of mesh points
        INTEGER, DIMENSION(:),     INTENT(IN   ) :: coord
        INTEGER,                   INTENT(IN   ) :: oldlabel
        INTEGER,                   INTENT(IN   ) :: newlabel
        INTEGER,                   INTENT(IN   ) :: Sm
        !!! mesh starting point number whether it will start from 0(ghost) or 1
        INTEGER,                   INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_rc_list) :: seedlst
        TYPE(ppm_rc_list) :: seedlsti

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(:),   POINTER     :: seedn
        INTEGER, DIMENSION(__DIME)                :: ld
#if   __DIME == __3D
        INTEGER, DIMENSION(__DIME)                :: ldi
#endif
        INTEGER, DIMENSION(__DIME)                :: Nm,Nmm
        INTEGER                              :: xx,xs,xe
        INTEGER                              :: Smp

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_floodFillScanline"

        LOGICAL :: south,north
#if   __DIME == __3D
        LOGICAL :: top,bottom
#endif
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (newlabel.EQ.oldlabel) GOTO 9999

        Nm=Nm_(1:__DIME)+MERGE(1,0,Sm.EQ.0)
        Nmm=Nm-1

        Smp=Sm+1

#if   __DIME == __2D
        IF (ABS(labels_(coord(1),coord(2))).EQ.newlabel) THEN
           labels_(coord(1),coord(2))=SIGN(oldlabel,labels_(coord(1),coord(2)))
#elif __DIME == __3D
        IF (ABS(labels_(coord(1),coord(2),coord(3))).EQ.newlabel) THEN
           labels_(coord(1),coord(2),coord(3))=SIGN(oldlabel,labels_(coord(1),coord(2),coord(3)))
#endif
        ENDIF

        !TODO
        !it only has been implemented for 6 and 4 connectivity forground
        !connectivity(3DIM and 2DIM respectively.)
        CALL seedlst%add(coord)
        SELECT CASE (FG_ConnectivityType%ID)
        CASE (1,3)
           seedlst_loop: DO WHILE (ASSOCIATED(seedlst%first))
              seedn => seedlst%first%getValue()
              ld=seedn
              CALL seedlst%remove()

              IF (ANY(ld.LT.Sm.OR.ld.GT.Nm)) THEN
                 CYCLE seedlst_loop
              ENDIF

#if   __DIME == __2D
              DO xs=ld(1),Sm,-1
                 IF (ABS(labels_(xs,ld(2))).NE.oldlabel) EXIT
              ENDDO
              xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.ABS(labels_(Sm,ld(2))).EQ.oldlabel)

              DO xe=ld(1)+1,Nm(1)
                 IF (ABS(labels_(xe,ld(2))).NE.oldlabel) EXIT
              ENDDO
              xe=xe-1

              south=.FALSE.
              north=.FALSE.
              IF      (ld(2).EQ.Sm) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))

                    IF (.NOT.north.AND.ABS(labels_(xx,Smp)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,Smp)
                       north=.TRUE.
                    ELSE IF (north.AND.ABS(labels_(xx,Smp)).NE.oldlabel) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE IF (ld(2).EQ.Nm(2)) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))

                    IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2))).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,Nmm(2))
                       south=.TRUE.
                    ELSE IF (south.AND.ABS(labels_(xx,Nmm(2))).NE.oldlabel) THEN
                       south=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE
                 DO xx=xs,xe
                    labels_(xx,ld(2))=SIGN(newlabel,labels_(xx,ld(2)))

                    IF (.NOT.south.AND.ABS(labels_(xx,ld(2)-1)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,ld(2)-1)
                       south=.TRUE.
                    ELSE IF (south.AND.ABS(labels_(xx,ld(2)-1)).NE.oldlabel) THEN
                       south=.FALSE.
                    ENDIF
                    IF (.NOT.north.AND.ABS(labels_(xx,ld(2)+1)).EQ.oldlabel) THEN
                       CALL seedlst%add(xx,ld(2)+1)
                       north=.TRUE.
                    ELSE IF (north.AND.ABS(labels_(xx,ld(2)+1)).NE.oldlabel) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ENDIF
#elif __DIME == __3D
              CALL seedlsti%add(ld)
              seedlsti_loop: DO WHILE (ASSOCIATED(seedlsti%first))
                 seedn => seedlsti%first%getValue()
                 ldi=seedn
                 CALL seedlsti%remove()

                 IF (ANY(ldi.LT.Sm.OR.ldi.GT.Nm)) THEN
                    CYCLE seedlsti_loop
                 ENDIF

                 DO xs=ldi(1),Sm,-1
                    IF (ABS(labels_(xs,ldi(2),ldi(3))).NE.oldlabel) EXIT
                 ENDDO
                 xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.ABS(labels_(Sm,ldi(2),ldi(3))).EQ.oldlabel)

                 DO xe=ldi(1)+1,Nm(1)
                    IF (ABS(labels_(xe,ldi(2),ldi(3))).NE.oldlabel) EXIT
                 ENDDO
                 xe=xe-1

                 south =.FALSE.
                 north =.FALSE.
                 bottom=.FALSE.
                 top   =.FALSE.
                 IF      (ldi(3).EQ.Sm) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,Sm)).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,Sm)
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,Sm)).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,Sm,Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,Sm,Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,Sm,Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),Smp)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),Smp)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !(ldi(2).EQ.Sm)
                 ELSE IF (ldi(3).EQ.Nm(3)) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),Nmm(3))).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ELSE
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.north.AND.ABS(labels_(xx,Smp,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,Smp,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,Nmm(2),ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=SIGN(newlabel,labels_(xx,ldi(2),ldi(3)))

                          IF (.NOT.south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.ABS(labels_(xx,ldi(2)-1,ldi(3))).NE.oldlabel) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).EQ.oldlabel) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.ABS(labels_(xx,ldi(2)+1,ldi(3))).NE.oldlabel) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.ABS(labels_(xx,ldi(2),ldi(3)-1)).NE.oldlabel) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).EQ.oldlabel) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.ABS(labels_(xx,ldi(2),ldi(3)+1)).NE.oldlabel) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ENDIF !ldi(3).EQ.Sm
              ENDDO seedlsti_loop !ASSOCIATED(seedlsti%first)
#endif
           ENDDO seedlst_loop !ASSOCIATED(seedlst%first)

        CASE DEFAULT
           fail("It has not implemented for this Foreground connectivity",ppm_error=ppm_error_fatal)

        END SELECT

        CALL seedlst%destroy()
        CALL seedlsti%destroy()

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_floodFillScanline_)


      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_floodfill
      !-------------------------------------------------------------------------
      !
      !  Purpose      : The purpose of floodFill is to color an entire area
      !                 of connected pixels with the same color
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       :
      !                 info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Author           - y.afshar           April   2014
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_floodFillScanlineConditional)(labels_,  &
      &          field_,Nm_,coord,oldlabel,newlabel,lower_,upper_,Sm,info)
        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: labels_
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
#endif

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: field_
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: field_
#endif

        INTEGER, DIMENSION(:),     POINTER       :: Nm_
        !!! number of mesh points
        INTEGER, DIMENSION(:),     INTENT(IN   ) :: coord
        INTEGER,                   INTENT(IN   ) :: oldlabel
        INTEGER,                   INTENT(IN   ) :: newlabel

        REAL(MK),                  INTENT(IN   ) :: lower_
        REAL(MK),                  INTENT(IN   ) :: upper_

        INTEGER,                   INTENT(IN   ) :: Sm
        !!! mesh starting point number whether it will start from 0(ghost) or 1
        INTEGER,                   INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_rc_list) :: seedlst
        TYPE(ppm_rc_list) :: seedlsti

        REAL(ppm_kind_double) :: t0
        REAL(MK) :: lower,upper

        INTEGER, DIMENSION(:),   POINTER :: seedn
        INTEGER, DIMENSION(__DIME)                :: ld
#if   __DIME == __3D
        INTEGER, DIMENSION(__DIME)                :: ldi
#endif
        INTEGER, DIMENSION(__DIME)            :: Nm,Nmm
        INTEGER                          :: xx,xs,xe
        INTEGER                          :: Smp

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_floodFillScanline"

        LOGICAL :: south,north
#if   __DIME == __3D
        LOGICAL :: top,bottom
#endif
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (newlabel.EQ.oldlabel) GOTO 9999

        Nm=Nm_(1:__DIME)+MERGE(1,0,Sm.EQ.0)
        Nmm=Nm-1

        Smp=Sm+1

#if   __DIME == __2D
        IF (ABS(labels_(coord(1),coord(2))).EQ.newlabel) THEN
           labels_(coord(1),coord(2))=SIGN(oldlabel,labels_(coord(1),coord(2)))
#elif __DIME == __3D
        IF (ABS(labels_(coord(1),coord(2),coord(3))).EQ.newlabel) THEN
           labels_(coord(1),coord(2),coord(3))=SIGN(oldlabel,labels_(coord(1),coord(2),coord(3)))
#endif
        ENDIF

        CALL seedlst%add(coord)
        SELECT CASE (FG_ConnectivityType%ID)
        CASE (1,3)
           seedlst_loop: DO WHILE (ASSOCIATED(seedlst%first))
              seedn => seedlst%first%getValue()
              ld=seedn
              CALL seedlst%remove()

              IF (ANY(ld.LT.Sm.OR.ld.GT.Nm)) THEN
                 CYCLE seedlst_loop
              ENDIF

#if   __DIME == __2D
              lower=field_(ld(1),ld(2))-lower_
              upper=field_(ld(1),ld(2))+upper_

              DO xs=ld(1),Sm,-1
                 IF (labels_(xs,ld(2)).NE.oldlabel.OR. &
                 &   field_(xs,ld(2)).LT.lower.OR.field_(xs,ld(2)).GT.upper) EXIT
              ENDDO
              xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.labels_(Sm,ld(2)).EQ.oldlabel &
              & .AND.lower.LE.field_(Sm,ld(2)).AND.field_(Sm,ld(2)).LE.upper)

              DO xe=ld(1)+1,Nm(1)
                 IF (labels_(xe,ld(2)).NE.oldlabel.OR. &
                 & field_(xe,ld(2)).LT.lower.OR.field_(xe,ld(2)).GT.upper) EXIT
              ENDDO
              xe=xe-1

              south=.FALSE.
              north=.FALSE.
              IF      (ld(2).EQ.Sm) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=newlabel

                    IF (.NOT.north.AND.labels_(xx,Smp).EQ.oldlabel &
                    &   .AND.lower.LE.field_(xx,Smp).AND.field_(xx,Smp).LE.upper) THEN
                       CALL seedlst%add(xx,Smp)
                       north=.TRUE.
                    ELSE IF (north.AND.labels_(xx,Smp).NE.oldlabel &
                    &   .AND.(field_(xx,Smp).LT.lower.OR.field_(xx,Smp).GT.upper)) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE IF (ld(2).EQ.Nm(2)) THEN
                 DO xx=xs,xe
                    labels_(xx,ld(2))=newlabel

                    IF (.NOT.south.AND.labels_(xx,Nmm(2)).EQ.oldlabel &
                    &   .AND.lower.LE.field_(xx,Nmm(2)).AND.field_(xx,Nmm(2)).LE.upper) THEN
                       CALL seedlst%add(xx,Nmm(2))
                       south=.TRUE.
                    ELSE IF (south.AND.labels_(xx,Nmm(2)).NE.oldlabel &
                    &   .AND.(field_(xx,Nmm(2)).LT.lower.OR.field_(xx,Nmm(2)).GT.upper)) THEN
                       south=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ELSE
                 DO xx=xs,xe
                    labels_(xx,ld(2))=newlabel

                    IF (.NOT.south.AND.labels_(xx,ld(2)-1).EQ.oldlabel &
                    &   .AND.lower.LE.field_(xx,ld(2)-1).AND.field_(xx,ld(2)-1).LE.upper) THEN
                       CALL seedlst%add(xx,ld(2)-1)
                       south=.TRUE.
                    ELSE IF (south.AND.labels_(xx,ld(2)-1).NE.oldlabel &
                    &   .AND.(field_(xx,ld(2)-1).LT.lower.OR.field_(xx,ld(2)-1).GT.upper)) THEN
                       south=.FALSE.
                    ENDIF
                    IF (.NOT.north.AND.labels_(xx,ld(2)+1).EQ.oldlabel &
                    &   .AND.lower.LE.field_(xx,ld(2)+1).AND.field_(xx,ld(2)+1).LE.upper) THEN
                       CALL seedlst%add(xx,ld(2)+1)
                       north=.TRUE.
                    ELSE IF (north.AND.labels_(xx,ld(2)+1).NE.oldlabel &
                    &   .AND.(field_(xx,ld(2)+1).LT.lower.OR.field_(xx,ld(2)+1).GT.upper)) THEN
                       north=.FALSE.
                    ENDIF
                 ENDDO !xx=xs,xe
              ENDIF
#elif __DIME == __3D
              CALL seedlsti%add(ld)
              seedlsti_loop: DO WHILE (ASSOCIATED(seedlsti%first))
                 seedn => seedlsti%first%getValue()
                 ldi=seedn
                 CALL seedlsti%remove()

                 IF (ANY(ldi.LT.Sm.OR.ldi.GT.Nm)) THEN
                    CYCLE seedlsti_loop
                 ENDIF

                 lower=field_(ldi(1),ldi(2),ldi(3))-lower_
                 upper=field_(ldi(1),ldi(2),ldi(3))+upper_

                 DO xs=ldi(1),Sm,-1
                    IF (labels_(xs,ldi(2),ldi(3)).NE.oldlabel.OR. &
                    &  field_(xs,ldi(2),ldi(3)).LT.lower.OR.field_(xs,ldi(2),ldi(3)).GT.upper) EXIT
                 ENDDO
                 xs=MERGE(Sm,xs+1,xs.EQ.Sm.AND.labels_(Sm,ldi(2),ldi(3)).EQ.oldlabel.AND. &
                 & lower.LE.field_(Sm,ldi(2),ldi(3)).AND.field_(Sm,ldi(2),ldi(3)).LE.upper)

                 DO xe=ldi(1)+1,Nm(1)
                    IF (labels_(xe,ldi(2),ldi(3)).NE.oldlabel.OR. &
                    &  field_(xe,ldi(2),ldi(3)).LT.lower.OR.field_(xe,ldi(2),ldi(3)).GT.upper) EXIT
                 ENDDO
                 xe=xe-1

                 south =.FALSE.
                 north =.FALSE.
                 bottom=.FALSE.
                 top   =.FALSE.
                 IF      (ldi(3).EQ.Sm) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.north.AND.labels_(xx,Smp,Sm).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Smp,Sm) &
                          &   .AND.field_(xx,Smp,Sm).LE.upper) THEN
                             CALL seedlsti%add(xx,Smp,Sm)
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,Smp,Sm).NE.oldlabel &
                          &   .AND.(field_(xx,Smp,Sm).LT.lower &
                          &   .OR.field_(xx,Smp,Sm).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.labels_(xx,Sm,Smp).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Sm,Smp) &
                          &   .AND.field_(xx,Sm,Smp).LE.upper) THEN
                             CALL seedlst%add(xx,Sm,Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,Sm,Smp).NE.oldlabel &
                          &   .AND.(field_(xx,Sm,Smp).LT.lower &
                          &   .OR.field_(xx,Sm,Smp).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,Nmm(2),ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Nmm(2),ldi(3)) &
                          &   .AND.field_(xx,Nmm(2),ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,Nmm(2),ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,Nmm(2),ldi(3)).LT.lower &
                          &   .OR.field_(xx,Nmm(2),ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.labels_(xx,ldi(2),Smp).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),Smp) &
                          &   .AND.field_(xx,ldi(2),Smp).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,ldi(2),Smp).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),Smp).LT.lower &
                          &   .OR.field_(xx,ldi(2),Smp).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,ldi(2)-1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)-1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)-1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,ldi(2)-1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)-1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)-1,ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.labels_(xx,ldi(2)+1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)+1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)+1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,ldi(2)+1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)+1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)+1,ldi(3)).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.top.AND.labels_(xx,ldi(2),Smp).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),Smp) &
                          &   .AND.field_(xx,ldi(2),Smp).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),Smp)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,ldi(2),Smp).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),Smp).LT.lower &
                          &   .OR.field_(xx,ldi(2),Smp).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !(ldi(2).EQ.Sm)
                 ELSE IF (ldi(3).EQ.Nm(3)) THEN
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.north.AND.labels_(xx,Smp,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Smp,ldi(3)) &
                          &   .AND.field_(xx,Smp,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,Smp,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,Smp,ldi(3)).LT.lower &
                          &   .OR.field_(xx,Smp,ldi(3)).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),Nmm(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),Nmm(3)) &
                          &   .AND.field_(xx,ldi(2),Nmm(3)).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),Nmm(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),Nmm(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2),Nmm(3)).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,Nmm(2),ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Nmm(2),ldi(3)) &
                          &   .AND.field_(xx,Nmm(2),ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,Nmm(2),ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,Nmm(2),ldi(3)).LT.lower &
                          &   .OR.field_(xx,Nmm(2),ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),Nmm(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),Nmm(3)) &
                          &   .AND.field_(xx,ldi(2),Nmm(3)).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),Nmm(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),Nmm(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2),Nmm(3)).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,ldi(2)-1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)-1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)-1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,ldi(2)-1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)-1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)-1,ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.labels_(xx,ldi(2)+1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)+1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)+1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,ldi(2)+1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)+1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)+1,ldi(3)).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),Nmm(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),Nmm(3)) &
                          &   .AND.field_(xx,ldi(2),Nmm(3)).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),Nmm(3))
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),Nmm(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),Nmm(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2),Nmm(3)).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ELSE
                    IF      (ldi(2).EQ.Sm) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.north.AND.labels_(xx,Smp,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Smp,ldi(3)) &
                          &   .AND.field_(xx,Smp,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,Smp,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,Smp,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,Smp,ldi(3)).LT.lower &
                          &   .OR.field_(xx,Smp,ldi(3)).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),ldi(3)-1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)-1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)-1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),ldi(3)-1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)-1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)-1).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.labels_(xx,ldi(2),ldi(3)+1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)+1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)+1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,ldi(2),ldi(3)+1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)+1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)+1).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE IF (ldi(2).EQ.Nm(2)) THEN
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,Nmm(2),ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,Nmm(2),ldi(3)) &
                          &   .AND.field_(xx,Nmm(2),ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,Nmm(2),ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,Nmm(2),ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,Nmm(2),ldi(3)).LT.lower &
                          &   .OR.field_(xx,Nmm(2),ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),ldi(3)-1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)-1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)-1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),ldi(3)-1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)-1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)-1).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.labels_(xx,ldi(2),ldi(3)+1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)+1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)+1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,ldi(2),ldi(3)+1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)+1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)+1).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ELSE
                       DO xx=xs,xe
                          labels_(xx,ldi(2),ldi(3))=newlabel

                          IF (.NOT.south.AND.labels_(xx,ldi(2)-1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)-1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)-1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)-1,ldi(3))
                             south=.TRUE.
                          ELSE IF (south.AND.labels_(xx,ldi(2)-1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)-1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)-1,ldi(3)).GT.upper)) THEN
                             south=.FALSE.
                          ENDIF
                          IF (.NOT.north.AND.labels_(xx,ldi(2)+1,ldi(3)).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2)+1,ldi(3)) &
                          &   .AND.field_(xx,ldi(2)+1,ldi(3)).LE.upper) THEN
                             CALL seedlsti%add(xx,ldi(2)+1,ldi(3))
                             north=.TRUE.
                          ELSE IF (north.AND.labels_(xx,ldi(2)+1,ldi(3)).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2)+1,ldi(3)).LT.lower &
                          &   .OR.field_(xx,ldi(2)+1,ldi(3)).GT.upper)) THEN
                             north=.FALSE.
                          ENDIF

                          IF (.NOT.bottom.AND.labels_(xx,ldi(2),ldi(3)-1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)-1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)-1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)-1)
                             bottom=.TRUE.
                          ELSE IF (bottom.AND.labels_(xx,ldi(2),ldi(3)-1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)-1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)-1).GT.upper)) THEN
                             bottom=.FALSE.
                          ENDIF
                          IF (.NOT.top.AND.labels_(xx,ldi(2),ldi(3)+1).EQ.oldlabel &
                          &   .AND.lower.LE.field_(xx,ldi(2),ldi(3)+1) &
                          &   .AND.field_(xx,ldi(2),ldi(3)+1).LE.upper) THEN
                             CALL seedlst%add(xx,ldi(2),ldi(3)+1)
                             top=.TRUE.
                          ELSE IF (top.AND.labels_(xx,ldi(2),ldi(3)+1).NE.oldlabel &
                          &   .AND.(field_(xx,ldi(2),ldi(3)+1).LT.lower &
                          &   .OR.field_(xx,ldi(2),ldi(3)+1).GT.upper)) THEN
                             top=.FALSE.
                          ENDIF
                       ENDDO !xx=xs,xe
                    ENDIF !ldi(2).EQ.Sm
                 ENDIF !ldi(3).EQ.Sm
              ENDDO seedlsti_loop !ASSOCIATED(seedlsti%first)
#endif
           ENDDO seedlst_loop !ASSOCIATED(seedlst%first)

        END SELECT

        CALL seedlst%destroy()
        CALL seedlsti%destroy()

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_floodFillScanlineConditional)
