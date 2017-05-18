#if   __DIME == __3D
        SUBROUTINE E_PS_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PS)            :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("E_PS_destroy")

          DEALLOCATE(this%NeighborsPoints,STAT=info)
          or_fail_dealloc("this%NeighborsPoints")

          end_subroutine()

        END SUBROUTINE E_PS_destroy
#endif

        SUBROUTINE DTYPE(E_PS_PrepareEnergyCalculation)(this,info)
          !!! Method is used to prepare energy functions.
          !!! It is called only once in the beginning of the filter.

          IMPLICIT NONE

          CLASS(E_PS)            :: this

          INTEGER, INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: s,sX,sY
#if   __DIME == __3D
          REAL(MK) :: sZ
#endif

          INTEGER :: sizeX,sizeY
          INTEGER :: i,j,lX,lY,l
#if   __DIME == __3D
          INTEGER :: sizeZ
          INTEGER :: k,lZ
#endif

          start_subroutine("PrepareEnergyCalculation")

          IF (ALLOCATED(this%NeighborsPoints)) THEN
             DEALLOCATE(this%NeighborsPoints,STAT=info)
             or_fail_dealloc("this%NeighborsPoints")
          ENDIF

          !!! the radius is expected to be given in px size of the first
          !!! axis. We scale for the all the following dimensions according
          !!! to the image spacing.
          e_dX=CEILING(this%m_Radius)
          !!! scaling factor
          sX=REAL(e_dX*e_dX,MK)

          e_dY=CEILING(this%m_Radius*pixel(1)/pixel(2))
          sY=REAL(e_dY*e_dY,MK)

          !!! this index is used for index shift in the MASK array
          lX=e_dX+1
          lY=e_dY+1

          !!! the size of the MASK array
          sizeX=e_dX*2+1
          sizeY=e_dY*2+1

          l=0

#if   __DIME == __2D
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                ENDIF
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                   this%NeighborsPoints(1,l)=i+lX
                   this%NeighborsPoints(2,l)=j+lY
                ENDIF
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#elif __DIME == __3D
          e_dZ=CEILING(this%m_Radius*pixel(1)/pixel(3))
          sZ=REAL(e_dZ*e_dZ,MK)

          !!! this index is used for index shift in the MASK array
          lZ=e_dZ+1

          !!! the size of the MASK array
          sizeZ=e_dZ*2+1

          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+REAL(j*j)/sY+REAL(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+(j*j)/sY+(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                      this%NeighborsPoints(1,l)=i+lX
                      this%NeighborsPoints(2,l)=j+lY
                      this%NeighborsPoints(3,l)=k+lZ
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#endif

          end_subroutine()

        END SUBROUTINE DTYPE(E_PS_PrepareEnergyCalculation)

        FUNCTION DTYPE(E_PS_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,  &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PS)                                           :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PS_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d1,d2
          REAL(MK)                            :: intensity
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(ppm_kind_double)               :: vMeanFrom,vMeanTo

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom,vNTo

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PS_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom=-intensity
          vSumTo  =zero

          vNFrom=-1
          vNTo  =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom=vSumFrom + tmpimage(i,j)
                vNFrom  =vNFrom   + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo  =vSumTo   + tmpimage(i,j)
                vNTo    =vNTo     + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom=vSumFrom + tmpimage(i,j,k)
                vNFrom  =vNFrom   + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo  =vSumTo   + tmpimage(i,j,k)
                vNTo    =vNTo     + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo.EQ.0) THEN !this should only happen with the BG label
             vMeanTo=this%gSums(newlabel_region)/this%gCount(newlabel_region)
             !this%gCount(newlabel_region) can not be less than small
          ELSE
             vMeanTo=REAL(vSumTo,ppm_kind_double)/REAL(vNTo,ppm_kind_double)
          ENDIF

          IF (vNFrom.EQ.0) THEN
             vMeanFrom=this%gSums(oldlabel_region)/this%gCount(oldlabel_region)
             !this%gCount(oldlabel_region) can not be less than or equal to one
          ELSE
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/REAL(vNFrom,ppm_kind_double)
          ENDIF

          d1=REAL(intensity,ppm_kind_double)-vMeanTo
          d2=REAL(intensity,ppm_kind_double)-vMeanFrom
          d1=d1*d1
          d2=d2*d2

          DTYPE(E_PS_EvaluateEnergyDifference)=REAL(this%m_Coefficient*(d1-d2),MK)

        END FUNCTION DTYPE(E_PS_EvaluateEnergyDifference)

        FUNCTION DTYPE(E_PS_CalculateTotalEnergy)(this,imageIn,labelsIn,Nm,info)

          IMPLICIT NONE

          CLASS(E_PS)                                           :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: Nm
          INTEGER,                                INTENT(  OUT) :: info

          REAL(MK)                                              :: DTYPE(E_PS_CalculateTotalEnergy)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d1,tot,Mean
          REAL(ppm_kind_double)               :: intensity
          REAL(MK)                            :: vSum

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k,kk,lZ
#endif
          INTEGER                            :: i,j,l,ii,jj,lX,lY
          INTEGER                            :: vN,vlabel
          INTEGER, DIMENSION(__DIME)         :: istart
          INTEGER, DIMENSION(__DIME)         :: iend

          tot=zerod

          istart(1:__DIME) = Nm(__DIME+1:2*__DIME)
          iend(1:__DIME)   = Nm(2*__DIME+1:3*__DIME)

          lX=-e_dX-1
          lY=-e_dY-1
#if   __DIME == __2D
          DO j=1,Nm(2)
             DO i=1,Nm(1)
                intensity=REAL(imageIn(i,j),ppm_kind_double)
                tmpimage  =>  imageIn(i-e_dX:i+e_dX,j-e_dY:j+e_dY)
                vlabel=ABS(labelsIn(i,j))
                tmplabels => labelsIn(i-e_dX:i+e_dX,j-e_dY:j+e_dY)

                vSum=zero
                vN  =0

                IF (vlabel.EQ.0.OR.vlabel.EQ.FORBIDDEN) THEN
                   DO l=1,this%NumberOfNeighbors
                      ii=this%NeighborsPoints(1,l)
                      jj=this%NeighborsPoints(2,l)
                      IF (istart(1).EQ.1) THEN
                         IF (i+ii+lX.LT.1) CYCLE
                      ENDIF
                      IF (iend(1).EQ.Ngrid(1)) THEN
                         IF (i+ii+lX.GT.Nm(1)) CYCLE
                      ENDIF
                      IF (istart(2).EQ.1) THEN
                         IF (j+jj+lY.LT.1) CYCLE
                      ENDIF
                      IF (iend(2).EQ.Ngrid(2)) THEN
                         IF (j+jj+lY.GT.Nm(2)) CYCLE
                      ENDIF
                      vlabel=ABS(tmplabels(ii,jj))
                      IF (vlabel.EQ.0.OR.vlabel.EQ.FORBIDDEN) THEN
                         vSum=vSum + tmpimage(ii,jj)
                         vN  =vN   + 1
                      ENDIF
                   ENDDO
                ELSE
                   DO l=1,this%NumberOfNeighbors
                      ii=this%NeighborsPoints(1,l)
                      jj=this%NeighborsPoints(2,l)
                      IF (istart(1).EQ.1) THEN
                         IF (i+ii+lX.LT.1) CYCLE
                      ENDIF
                      IF (iend(1).EQ.Ngrid(1)) THEN
                         IF (i+ii+lX.GT.Nm(1)) CYCLE
                      ENDIF
                      IF (istart(2).EQ.1) THEN
                         IF (j+jj+lY.LT.1) CYCLE
                      ENDIF
                      IF (iend(2).EQ.Ngrid(2)) THEN
                         IF (j+jj+lY.GT.Nm(2)) CYCLE
                      ENDIF
                      IF (ABS(tmplabels(ii,jj)).EQ.vlabel) THEN
                         vSum=vSum + tmpimage(ii,jj)
                         vN  =vN   + 1
                      ENDIF
                   ENDDO
                ENDIF

                Mean=REAL(vSum,ppm_kind_double)/REAL(vN,ppm_kind_double)

                d1=Mean-intensity

                tot=tot+d1*d1
             ENDDO
          ENDDO
#elif __DIME == __3D
          lZ=-e_dZ-1
          DO k=1,Nm(3)
             DO j=1,Nm(2)
                DO i=1,Nm(1)
                   intensity=REAL(imageIn(i,j,k),ppm_kind_double)

                   tmpimage  =>  imageIn(i-e_dX:i+e_dX,j-e_dY:j+e_dY,k-e_dZ:k+e_dZ)

                   tmplabels => labelsIn(i-e_dX:i+e_dX,j-e_dY:j+e_dY,k-e_dZ:k+e_dZ)

                   vlabel=ABS(labelsIn(i,j,k))

                   vSum=zero
                   vN  =0

                   IF (vlabel.EQ.0.OR.vlabel.EQ.FORBIDDEN) THEN
                      DO l=1,this%NumberOfNeighbors
                         ii=this%NeighborsPoints(1,l)
                         jj=this%NeighborsPoints(2,l)
                         kk=this%NeighborsPoints(3,l)

                         IF (istart(1).EQ.1) THEN
                            IF (i+ii+lX.LT.1) CYCLE
                         ENDIF
                         IF (iend(1).EQ.Ngrid(1)) THEN
                            IF (i+ii+lX.GT.Nm(1)) CYCLE
                         ENDIF
                         IF (istart(2).EQ.1) THEN
                            IF (j+jj+lY.LT.1) CYCLE
                         ENDIF
                         IF (iend(2).EQ.Ngrid(2)) THEN
                            IF (j+jj+lY.GT.Nm(2)) CYCLE
                         ENDIF
                         IF (istart(3).EQ.1) THEN
                            IF (k+kk+lZ.LT.1) CYCLE
                         ENDIF
                         IF (iend(3).EQ.Ngrid(3)) THEN
                            IF (k+kk+lZ.GT.Nm(3)) CYCLE
                         ENDIF

                         vlabel=ABS(tmplabels(ii,jj,kk))
                         IF (vlabel.EQ.0.OR.vlabel.EQ.FORBIDDEN) THEN
                            vSum=vSum + tmpimage(ii,jj,kk)
                            vN  =vN   + 1
                         ENDIF
                      ENDDO
                   ELSE
                      DO l=1,this%NumberOfNeighbors
                         ii=this%NeighborsPoints(1,l)
                         jj=this%NeighborsPoints(2,l)
                         kk=this%NeighborsPoints(3,l)

                         IF (istart(1).EQ.1) THEN
                            IF (i+ii+lX.LT.1) CYCLE
                         ENDIF
                         IF (iend(1).EQ.Ngrid(1)) THEN
                            IF (i+ii+lX.GT.Nm(1)) CYCLE
                         ENDIF
                         IF (istart(2).EQ.1) THEN
                            IF (j+jj+lY.LT.1) CYCLE
                         ENDIF
                         IF (iend(2).EQ.Ngrid(2)) THEN
                            IF (j+jj+lY.GT.Nm(2)) CYCLE
                         ENDIF
                         IF (istart(3).EQ.1) THEN
                            IF (k+kk+lZ.LT.1) CYCLE
                         ENDIF
                         IF (iend(3).EQ.Ngrid(3)) THEN
                            IF (k+kk+lZ.GT.Nm(3)) CYCLE
                         ENDIF

                         IF (ABS(tmplabels(ii,jj,kk)).EQ.vlabel) THEN
                            vSum=vSum + tmpimage(ii,jj,kk)
                            vN  =vN   + 1
                         ENDIF
                      ENDDO
                   ENDIF

                   Mean=REAL(vSum,ppm_kind_double)/REAL(vN,ppm_kind_double)

                   d1=Mean-intensity

                   tot=tot+d1*d1
                ENDDO
             ENDDO
          ENDDO
#endif

          DTYPE(E_PS_CalculateTotalEnergy)=REAL(tot,MK)
          info=0

        END FUNCTION DTYPE(E_PS_CalculateTotalEnergy)

        FUNCTION DTYPE(E_PS_EvaluateEnergyDifference_E_Merge)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,          &
        &        oldlabel_region,newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(E_PS)                                           :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(E_PS_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d1,d2
          REAL(MK)                            :: intensity,tmp
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(ppm_kind_double)               :: vNTo,vNFrom
          REAL(MK)                            :: vSumOfSqFrom,vSumOfSqTo
          REAL(ppm_kind_double)               :: vMeanTo,vMeanFrom
          REAL(ppm_kind_double)               :: vVarTo,vVarFrom

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom_,vNTo_

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom    =-intensity
          vSumOfSqFrom=-intensity*intensity
          vSumTo      =zero
          vSumOfSqTo  =zero

          vNFrom_     =-1
          vNTo_       =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom     = vSumFrom     + tmpimage(i,j)
                vSumOfSqFrom = vSumOfSqFrom + tmpimage(i,j)*tmpimage(i,j)

                vNFrom_      = vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo     = vSumTo     + tmpimage(i,j)
                vSumOfSqTo = vSumOfSqTo + tmpimage(i,j)*tmpimage(i,j)

                vNTo_      = vNTo_      + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom     = vSumFrom     + tmpimage(i,j,k)
                vSumOfSqFrom = vSumOfSqFrom + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNFrom_      = vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo     = vSumTo     + tmpimage(i,j,k)
                vSumOfSqTo = vSumOfSqTo + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNTo_      = vNTo_      + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo_.EQ.0) THEN !this should only happen with the BG label
             !this%gCount(newlabel_region) can not be less than small
             vNTo=this%gCount(newlabel_region)
             vMeanTo=this%gSums(newlabel_region)/vNTo
             vVarTo=this%CalculateVariance(this%gSumsq(newlabel_region), &
             &      vMeanTo,vNTo)
          ELSE
             vNTo=REAL(vNTo_,ppm_kind_double)
             vMeanTo=REAL(vSumTo,ppm_kind_double)/vNTo
             vVarTo=this%CalculateVariance(REAL(vSumOfSqTo,ppm_kind_double),vMeanTo,vNTo)
          ENDIF

          IF (vNFrom_.EQ.0) THEN
             !this%gCount(oldlabel_region) can not be less than or equal one
             vNFrom=this%gCount(oldlabel_region)
             vMeanFrom=this%gSums(oldlabel_region)/vNFrom
             vVarFrom=this%CalculateVariance(this%gSumsq(oldlabel_region), &
             &        vMeanFrom,vNFrom)
          ELSE
             vNFrom=REAL(vNFrom_,ppm_kind_double)
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/vNFrom

             vVarFrom=this%CalculateVariance(REAL(vSumOfSqFrom,ppm_kind_double),vMeanFrom,vNFrom)
          ENDIF

          IF (oldlabel.NE.0.AND.newlabel.NE.0) THEN
             tmp=this%CalculateKullbackLeiblerDistance(vMeanTo, &
             &   vMeanFrom,vVarTo,vVarFrom,vNTo,vNFrom)
             e_merge=tmp.LT.this%RegionMergingThreshold
          ELSE
             e_merge=.FALSE.
          ENDIF

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PS_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ENDIF

          d1=REAL(intensity,ppm_kind_double)-vMeanTo
          d2=REAL(intensity,ppm_kind_double)-vMeanFrom
          d1=d1*d1
          d2=d2*d2

          DTYPE(E_PS_EvaluateEnergyDifference_E_Merge)=REAL(this%m_Coefficient*(d1-d2),MK)

        END FUNCTION DTYPE(E_PS_EvaluateEnergyDifference_E_Merge)

#if   __DIME == __3D
        SUBROUTINE E_PSGaussian_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PSGaussian)    :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("E_PSGaussian_destroy")

          DEALLOCATE(this%NeighborsPoints,STAT=info)
          or_fail_dealloc("this%NeighborsPoints")

          end_subroutine()

        END SUBROUTINE E_PSGaussian_destroy
#endif

        SUBROUTINE DTYPE(E_PSGaussian_PrepareEnergyCalculation)(this,info)
          !!! Method is used to prepare energy functions.
          !!! It is called only once in the beginning of the filter.

          IMPLICIT NONE

          CLASS(E_PSGaussian)    :: this

          INTEGER, INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: s,sX,sY
#if   __DIME == __3D
          REAL(MK) :: sZ
#endif

          INTEGER :: sizeX,sizeY
          INTEGER :: i,j,lX,lY,l
#if   __DIME == __3D
          INTEGER :: sizeZ
          INTEGER :: k,lZ
#endif

          start_subroutine("PrepareEnergyCalculation")

          !!! the radius is expected to be given in px size of the first
          !!! axis. We scale for the all the following dimensions according
          !!! to the image spacing.
          e_dX=CEILING(this%m_Radius)
          !!! scaling factor
          sX=REAL(e_dX*e_dX,MK)

          e_dY=CEILING(this%m_Radius*pixel(1)/pixel(2))
          sY=REAL(e_dY*e_dY,MK)

          !!! this index is used for index shift in the MASK array
          lX=e_dX+1
          lY=e_dY+1

          !!! the size of the MASK array
          sizeX=e_dX*2+1
          sizeY=e_dY*2+1

          l=0

#if   __DIME == __2D
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                ENDIF
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                   this%NeighborsPoints(1,l)=i+lX
                   this%NeighborsPoints(2,l)=j+lY
                ENDIF
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#elif __DIME == __3D
          e_dZ=CEILING(this%m_Radius*pixel(1)/pixel(3))
          sZ=REAL(e_dZ*e_dZ,MK)

          !!! this index is used for index shift in the MASK array
          lZ=e_dZ+1

          !!! the size of the MASK array
          sizeZ=e_dZ*2+1

          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+REAL(j*j)/sY+REAL(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+(j*j)/sY+(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                      this%NeighborsPoints(1,l)=i+lX
                      this%NeighborsPoints(2,l)=j+lY
                      this%NeighborsPoints(3,l)=k+lZ
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#endif

          end_subroutine()

        END SUBROUTINE DTYPE(E_PSGaussian_PrepareEnergyCalculation)

        FUNCTION DTYPE(E_PSGaussian_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,          &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PSGaussian)                                   :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PSGaussian_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d1,d2
          REAL(MK)                            :: intensity
          REAL(ppm_kind_double)               :: vNTo,vNFrom
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(MK)                            :: vSumOfSqFrom,vSumOfSqTo
          REAL(ppm_kind_double)               :: vMeanTo,vMeanFrom
          REAL(ppm_kind_double)               :: vVarTo,vVarFrom

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom_,vNTo_

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom    =-intensity
          vSumOfSqFrom=-intensity*intensity
          vSumTo      =zero
          vSumOfSqTo  =zero

          vNFrom_     =-1
          vNTo_       =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j)*tmpimage(i,j)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j)*tmpimage(i,j)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j,k)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j,k)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo_.EQ.0) THEN !this should only happen with the BG label
             !this%gCount(newlabel_region) can not be less than small
             vNTo=this%gCount(newlabel_region)
             vMeanTo=this%gSums(newlabel_region)/vNTo
             vVarTo=this%CalculateVariance(this%gSumsq(newlabel_region), &
             &      vMeanTo,vNTo)
          ELSE
             vNTo=REAL(vNTo_,ppm_kind_double)
             vMeanTo=REAL(vSumTo,ppm_kind_double)/vNTo
             vVarTo=this%CalculateVariance(REAL(vSumOfSqTo,ppm_kind_double),vMeanTo,vNTo)
          ENDIF
          IF (vVarTo.LE.smalld) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

          IF (vNFrom_.EQ.0) THEN
             !this%gCount(oldlabel_region) can not be less than or equal one
             vNFrom=this%gCount(oldlabel_region)
             vMeanFrom=this%gSums(oldlabel_region)/vNFrom
             vVarFrom=this%CalculateVariance(this%gSumsq(oldlabel_region), &
             &        vMeanFrom,vNFrom)
          ELSE
             vNFrom=REAL(vNFrom_,ppm_kind_double)
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/vNFrom
             vVarFrom=this%CalculateVariance(REAL(vSumOfSqFrom,ppm_kind_double),vMeanFrom,vNFrom)
          ENDIF
          IF (vVarFrom.LE.smalld) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference)=bigs
             RETURN
          ENDIF

          d1=REAL(intensity,ppm_kind_double)-vMeanTo
          d2=REAL(intensity,ppm_kind_double)-vMeanFrom
          d1=d1*d1/(twod*vVarTo)  +halfd*LOG(twod*pid*vVarTo)
          d2=d2*d2/(twod*vVarFrom)+halfd*LOG(twod*pid*vVarFrom)

          DTYPE(E_PSGaussian_EvaluateEnergyDifference)=REAL(this%m_Coefficient*(d1-d2),MK)

        END FUNCTION DTYPE(E_PSGaussian_EvaluateEnergyDifference)


        FUNCTION DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,oldlabel_region,  &
        &        newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(E_PSGaussian)                                   :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,                DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d1,d2
          REAL(MK)                            :: intensity,tmp
          REAL(ppm_kind_double)               :: vNTo,vNFrom
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(MK)                            :: vSumOfSqFrom,vSumOfSqTo
          REAL(ppm_kind_double)               :: vMeanTo,vMeanFrom
          REAL(ppm_kind_double)               :: vVarTo,vVarFrom

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom_,vNTo_

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom    =-intensity
          vSumOfSqFrom=-intensity*intensity
          vSumTo      =zero
          vSumOfSqTo  =zero

          vNFrom_     =-1
          vNTo_       =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j)*tmpimage(i,j)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j)*tmpimage(i,j)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j,k)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j,k)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNTo_     = vNTo_     + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo_.EQ.0) THEN !this should only happen with the BG label
             !this%gCount(newlabel_region) can not be less than small
             vNTo=this%gCount(newlabel_region)
             vMeanTo=this%gSums(newlabel_region)/vNTo
             vVarTo=this%CalculateVariance(this%gSumsq(newlabel_region), &
             &      vMeanTo,vNTo)
          ELSE
             vNTo=REAL(vNTo_,ppm_kind_double)
             vMeanTo=REAL(vSumTo,ppm_kind_double)/vNTo
             vVarTo=this%CalculateVariance(REAL(vSumOfSqTo,ppm_kind_double),vMeanTo,vNTo)
          ENDIF


          IF (vNFrom_.EQ.0) THEN
             !this%gCount(oldlabel_region) can not be less than or equal to one
             vNFrom=this%gCount(oldlabel_region)
             vMeanFrom=this%gSums(oldlabel_region)/vNFrom
             vVarFrom=this%CalculateVariance(this%gSumsq(oldlabel_region), &
             &        vMeanFrom,vNFrom)
          ELSE
             vNFrom=REAL(vNFrom_,ppm_kind_double)
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/vNFrom
             vVarFrom=this%CalculateVariance(REAL(vSumOfSqFrom,ppm_kind_double),vMeanFrom,vNFrom)
          ENDIF

          IF (oldlabel.NE.0.AND.newlabel.NE.0) THEN
             tmp=this%CalculateKullbackLeiblerDistance(vMeanTo, &
             &   vMeanFrom,vVarTo,vVarFrom,vNTo,vNFrom)
             e_merge=tmp.LT.this%RegionMergingThreshold
          ELSE
             e_merge=.FALSE.
          ENDIF

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ELSE IF (vVarTo.LE.smalld) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ELSE IF (vVarFrom.LE.smalld) THEN
             DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)=bigs
             RETURN
          ENDIF

          d1=REAL(intensity,ppm_kind_double)-vMeanTo
          d2=REAL(intensity,ppm_kind_double)-vMeanFrom
          d1=d1*d1/(twod*vVarTo)  +halfd*LOG(twod*pid*vVarTo)
          d2=d2*d2/(twod*vVarFrom)+halfd*LOG(twod*pid*vVarFrom)

          DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)=REAL(this%m_Coefficient*(d1-d2),MK)

        END FUNCTION DTYPE(E_PSGaussian_EvaluateEnergyDifference_E_Merge)

#if   __DIME == __3D
        SUBROUTINE E_PSPoisson_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PSPoisson)    :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("E_PSPoisson_destroy")

          DEALLOCATE(this%NeighborsPoints,STAT=info)
          or_fail_dealloc("this%NeighborsPoints")

          end_subroutine()

        END SUBROUTINE E_PSPoisson_destroy
#endif

        SUBROUTINE DTYPE(E_PSPoisson_PrepareEnergyCalculation)(this,info)
          !!! Method is used to prepare energy functions.
          !!! It is called only once in the beginning of the filter.

          IMPLICIT NONE

          CLASS(E_PSPoisson)     :: this

          INTEGER, INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: s,sX,sY
#if   __DIME == __3D
          REAL(MK) :: sZ
#endif

          INTEGER :: sizeX,sizeY
          INTEGER :: i,j,lX,lY,l
#if   __DIME == __3D
          INTEGER :: sizeZ
          INTEGER :: k,lZ
#endif

          start_subroutine("PrepareEnergyCalculation")

          !!! the radius is expected to be given in px size of the first
          !!! axis. We scale for the all the following dimensions according
          !!! to the image spacing.
          e_dX=CEILING(this%m_Radius)
          !!! scaling factor
          sX=REAL(e_dX*e_dX,MK)

          e_dY=CEILING(this%m_Radius*pixel(1)/pixel(2))
          sY=REAL(e_dY*e_dY,MK)

          !!! this index is used for index shift in the MASK array
          lX=e_dX+1
          lY=e_dY+1

          !!! the size of the MASK array
          sizeX=e_dX*2+1
          sizeY=e_dY*2+1

          l=0

#if   __DIME == __2D
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                ENDIF
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO j=-e_dY,e_dY
             DO i=-e_dX,e_dX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                   this%NeighborsPoints(1,l)=i+lX
                   this%NeighborsPoints(2,l)=j+lY
                ENDIF
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#elif __DIME == __3D
          e_dZ=CEILING(this%m_Radius*pixel(1)/pixel(3))
          sZ=REAL(e_dZ*e_dZ,MK)

          !!! this index is used for index shift in the MASK array
          lZ=e_dZ+1

          !!! the size of the MASK array
          sizeZ=e_dZ*2+1

          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+REAL(j*j)/sY+REAL(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO k=-e_dZ,e_dZ
             DO j=-e_dY,e_dY
                DO i=-e_dX,e_dX
                   s=REAL(i*i)/sX+(j*j)/sY+(k*k)/sZ
                   IF (s.LE.oneplus) THEN
                      l=l+1
                      this%NeighborsPoints(1,l)=i+lX
                      this%NeighborsPoints(2,l)=j+lY
                      this%NeighborsPoints(3,l)=k+lZ
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !!! Filling the NeighborsPoints array based on the hyper sphere location
#endif

          end_subroutine()

        END SUBROUTINE DTYPE(E_PSPoisson_PrepareEnergyCalculation)

        FUNCTION DTYPE(E_PSPoisson_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,         &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PSPoisson)                                    :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,                DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PSPoisson_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d
          REAL(MK)                            :: intensity
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(MK)                            :: vSumOfSqFrom,vSumOfSqTo
          REAL(ppm_kind_double)               :: vNTo,vNFrom
          REAL(ppm_kind_double)               :: vMeanTo,vMeanFrom
          REAL(ppm_kind_double)               :: vVarTo,vVarFrom

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom_,vNTo_

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom    =-intensity
          vSumOfSqFrom=-intensity*intensity
          vSumTo      =zero
          vSumOfSqTo  =zero

          vNFrom_     =-1
          vNTo_       =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j)*tmpimage(i,j)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j)*tmpimage(i,j)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j,k)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j,k)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo_.EQ.0) THEN !this should only happen with the BG label
             !this%gCount(newlabel_region) can not be less than small
             vNTo=this%gCount(newlabel_region)
             vMeanTo=this%gSums(newlabel_region)/vNTo
             vVarTo=this%CalculateVariance(this%gSumsq(newlabel_region), &
             &      vMeanTo,vNTo)
          ELSE
             vNTo=REAL(vNTo_,ppm_kind_double)
             vMeanTo=REAL(vSumTo,ppm_kind_double)/vNTo
             vVarTo=this%CalculateVariance(REAL(vSumOfSqTo,ppm_kind_double),vMeanTo,vNTo)
          ENDIF


          IF (vNFrom_.EQ.0) THEN
             !this%gCount(oldlabel_region) can not be less than or equal to one
             vNFrom=this%gCount(oldlabel_region)
             vMeanFrom=this%gSums(oldlabel_region)/vNFrom
             vVarFrom=this%CalculateVariance(this%gSumsq(oldlabel_region), &
             &        vMeanFrom,vNFrom)
          ELSE
             vNFrom=REAL(vNFrom_,ppm_kind_double)
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/vNFrom
             vVarFrom=this%CalculateVariance(REAL(vSumOfSqFrom,ppm_kind_double),vMeanFrom,vNFrom)
          ENDIF

          IF (vMeanFrom.LT.smalld.AND.vMeanTo.LT.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference)=zero
             RETURN
          ELSE IF (vMeanFrom.LE.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference)=REAL(this%m_Coefficient*smallestd,MK)
             RETURN
          ELSE IF (vMeanTo.LE.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference)=bigs
             RETURN
          ENDIF

          d=REAL(intensity,ppm_kind_double)*LOG(vMeanFrom/vMeanTo)-vMeanFrom+vMeanTo

          DTYPE(E_PSPoisson_EvaluateEnergyDifference)=REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PSPoisson_EvaluateEnergyDifference)


        FUNCTION DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,oldlabel_region, &
        &        newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(E_PSPoisson)                                    :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
#endif
          INTEGER,                DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          REAL(MK), DIMENSION(:,:),   POINTER :: tmpimage
#elif __DIME == __3D
          REAL(MK), DIMENSION(:,:,:), POINTER :: tmpimage
#endif
          REAL(ppm_kind_double)               :: d
          REAL(MK)                            :: intensity,tmp
          REAL(MK)                            :: vSumFrom,vSumTo
          REAL(MK)                            :: vSumOfSqFrom,vSumOfSqTo
          REAL(ppm_kind_double)               :: vNTo,vNFrom
          REAL(ppm_kind_double)               :: vMeanTo,vMeanFrom
          REAL(ppm_kind_double)               :: vVarTo,vVarFrom

#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER                            :: k
#endif
          INTEGER                            :: i,j,l
          INTEGER                            :: vNFrom_,vNTo_

#if   __DIME == __2D
          intensity=imageIn(coord(1),coord(2))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY)
#elif __DIME == __3D
          intensity=imageIn(coord(1),coord(2),coord(3))

          tmpimage  =>  imageIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)

          tmplabels => labelsIn(coord(1)-e_dX:coord(1)+e_dX, &
          &                    coord(2)-e_dY:coord(2)+e_dY, &
          &                    coord(3)-e_dZ:coord(3)+e_dZ)
#endif

          vSumFrom    =-intensity
          vSumOfSqFrom=-intensity*intensity
          vSumTo      =zero
          vSumOfSqTo  =zero

          vNFrom_     =-1
          vNTo_       =0

#if   __DIME == __2D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j)*tmpimage(i,j)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j)*tmpimage(i,j)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#elif __DIME == __3D
          DO l=1,this%NumberOfNeighbors
             i=this%NeighborsPoints(1,l)
             j=this%NeighborsPoints(2,l)
             k=this%NeighborsPoints(3,l)
             IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                vSumFrom    =vSumFrom     + tmpimage(i,j,k)
                vSumOfSqFrom=vSumOfSqFrom + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNFrom_     =vNFrom_      + 1
             ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
                vSumTo    =vSumTo     + tmpimage(i,j,k)
                vSumOfSqTo=vSumOfSqTo + tmpimage(i,j,k)*tmpimage(i,j,k)

                vNTo_     =vNTo_      + 1
             ENDIF
          ENDDO
#endif

          IF (vNTo_.EQ.0) THEN !this should only happen with the BG label
             !this%gCount(newlabel_region) can not be less than small
             vNTo=this%gCount(newlabel_region)
             vMeanTo=this%gSums(newlabel_region)/vNTo
             vVarTo=this%CalculateVariance(this%gSumsq(newlabel_region), &
             &      vMeanTo,vNTo)
          ELSE
             vNTo=REAL(vNTo_,ppm_kind_double)
             vMeanTo=REAL(vSumTo,ppm_kind_double)/vNTo
             vVarTo=this%CalculateVariance(REAL(vSumOfSqTo,ppm_kind_double),vMeanTo,vNTo)
          ENDIF

          IF (vNFrom_.EQ.0) THEN
             !this%gCount(oldlabel_region) can not be less than or equal to one
             vNFrom=this%gCount(oldlabel_region)
             vMeanFrom=this%gSums(oldlabel_region)/vNFrom
             vVarFrom=this%CalculateVariance(this%gSumsq(oldlabel_region), &
             &        vMeanFrom,vNFrom)
          ELSE
             vNFrom=REAL(vNFrom_,ppm_kind_double)
             vMeanFrom=REAL(vSumFrom,ppm_kind_double)/vNFrom
             vVarFrom=this%CalculateVariance(REAL(vSumOfSqFrom,ppm_kind_double),vMeanFrom,vNFrom)
          ENDIF

          IF (oldlabel.NE.0.AND.newlabel.NE.0) THEN
             tmp=this%CalculateKullbackLeiblerDistance(vMeanTo, &
             &   vMeanFrom,vVarTo,vVarFrom,vNTo,vNFrom)
             e_merge=tmp.LT.this%RegionMergingThreshold
          ELSE
             e_merge=.FALSE.
          ENDIF

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ELSE IF (vMeanFrom.LT.smalld.AND.vMeanTo.LT.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ELSE IF (vMeanFrom.LE.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)=REAL(this%m_Coefficient*smallestd,MK)
             RETURN
          ELSE IF (vMeanTo.LE.smalld) THEN
             DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)=bigs
             RETURN
          ENDIF

          d=REAL(intensity,ppm_kind_double)*LOG(vMeanFrom/vMeanTo)-vMeanFrom+vMeanTo

          DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)=REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PSPoisson_EvaluateEnergyDifference_E_Merge)
