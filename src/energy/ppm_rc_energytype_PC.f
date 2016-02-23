        !!! Functions related to energy differences for a pixel change in a
        !!! piece-wise constant image model with different noise models.
        !!!
        !!! In case of 2 regions this image model is also called Chan-Vese model.
#if   __DIME == __3D
        SUBROUTINE E_PC_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PC)            :: this

          INTEGER, INTENT(  OUT) :: info

          info=0

        END SUBROUTINE E_PC_destroy
#endif

        FUNCTION DTYPE(E_PC_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,    &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PC)                                           :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER        :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER        :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PC_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vNewToMean,vOldFromMean,d1,d2
          REAL(ppm_kind_double) :: intensity

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PC_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          ! compute the new mean of the old region after this pixel has left
          vNewToMean  =this%gSums(newlabel_region) + intensity
          vNewToMean  =vNewToMean/(this%gCount(newlabel_region)+oned)
          vOldFromMean=this%gSums(oldlabel_region)/this%gCount(oldlabel_region)
          d1          =intensity-vNewToMean
          d2          =intensity-vOldFromMean
          d1          =d1*d1
          d2          =d2*d2

          DTYPE(E_PC_EvaluateEnergyDifference)=REAL(this%m_Coefficient*(d1-d2),MK)

        END FUNCTION DTYPE(E_PC_EvaluateEnergyDifference)

        FUNCTION DTYPE(E_PC_CalculateTotalEnergy)(this,image_,labels_,Nm,info)

          IMPLICIT NONE

          CLASS(E_PC)                                           :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: Nm
          INTEGER,                                INTENT(  OUT) :: info

          REAL(MK)                                              :: DTYPE(E_PC_CalculateTotalEnergy)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: tot,Mean,Mean0
          REAL(ppm_kind_double) :: d1,intensity

          INTEGER :: vlabel,label_region
          INTEGER :: i,j
#if   __DIME == __3D
          INTEGER :: k
#endif

          tot=zerod
          vlabel=0
          Mean0=this%gSums(0)/this%gCount(0)

#if   __DIME == __2D
          DO j=1,Nm(2)
             DO i=1,Nm(1)
                intensity=REAL(image_(i,j),ppm_kind_double)
                IF (labels_(i,j).EQ.0.OR.labels_(i,j).EQ.FORBIDDEN) THEN
                   d1=Mean0-intensity
                ELSE
                   IF (ABS(labels_(i,j)).NE.vlabel) THEN
                      vlabel=ABS(labels_(i,j))
                      label_region=htable%search(vlabel)
                      IF (label_region.EQ.htable_null) CYCLE
                      Mean=this%gSums(label_region)/this%gCount(label_region)
                   ENDIF
                   d1=Mean-intensity
                ENDIF
                tot=tot+d1*d1
             ENDDO
          ENDDO
#elif __DIME == __3D
          DO k=1,Nm(3)
             DO j=1,Nm(2)
                DO i=1,Nm(1)
                   intensity=REAL(image_(i,j,k),ppm_kind_double)
                   IF (labels_(i,j,k).EQ.0.OR.labels_(i,j,k).EQ.FORBIDDEN) THEN
                      d1=Mean0-intensity
                   ELSE
                      IF (ABS(labels_(i,j,k)).NE.vlabel) THEN
                         vlabel=ABS(labels_(i,j,k))
                         label_region=htable%search(vlabel)
                         IF (label_region.EQ.htable_null) CYCLE
                         Mean=this%gSums(label_region)/this%gCount(label_region)
                      ENDIF
                      d1=Mean-intensity
                   ENDIF
                   tot=tot+d1*d1
                ENDDO
             ENDDO
          ENDDO
#endif

          DTYPE(E_PC_CalculateTotalEnergy)=REAL(tot,MK)
          info=0

        END FUNCTION DTYPE(E_PC_CalculateTotalEnergy)

#if   __DIME == __3D
        SUBROUTINE E_PCGaussian_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PCGaussian)    :: this

          INTEGER, INTENT(  OUT) :: info

          info=0

        END SUBROUTINE E_PCGaussian_destroy
#endif

        FUNCTION DTYPE(E_PCGaussian_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,            &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PCGaussian)                                   :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PCGaussian_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double), PARAMETER :: vOneBySq2Pi = oned/SQRT(twod*pid)
          REAL(ppm_kind_double)            :: vNFrom,vNTo
          REAL(ppm_kind_double)            :: vNewFromMean,vOldFromMean
          REAL(ppm_kind_double)            :: vNewToMean,vOldToMean
          REAL(ppm_kind_double)            :: vNewFromVar,vOldFromVar
          REAL(ppm_kind_double)            :: vNewToVar,vOldToVar
          REAL(ppm_kind_double)            :: d
          REAL(ppm_kind_double)            :: intensity

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PCGaussian_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          vNFrom      =this%gCount(oldlabel_region)
          vNTo        =this%gCount(newlabel_region)

          vNewFromMean=(this%gSums(oldlabel_region)-intensity)/(vNFrom - oned)
          !vNFrom can not be less than or equal to one so it is safe to use this
          vOldFromMean=this%gSums(oldlabel_region)/vNFrom

          vNewToMean  =(this%gSums(newlabel_region)+intensity)/(vNTo + oned)
          vOldToMean  =this%gSums(newlabel_region)/vNTo
          !vNTo can not be zero, so it is safe to use this

          vNewFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region)- &
          &            intensity*intensity,vNewFromMean,vNFrom-oned)
          vOldFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region), &
          &            vOldFromMean,vNFrom)

          vNewToVar   =this%CalculateVariance(this%gSumsq(newlabel_region)+ &
          &            intensity*intensity,vNewToMean,vNTo+oned)
          vOldToVar   =this%CalculateVariance(this%gSumsq(newlabel_region), &
          &            vOldToMean,vNTo)

          ! it might happen that the variance is exactly 0 (if all pixel contain
          ! the same value). We therefore set it to a small value.
          IF (vNewFromVar.LT.smallestd) vNewFromVar=smalld*tend
          IF (vOldFromVar.LT.smallestd) vOldFromVar=smalld*tend
          IF (vNewToVar  .LT.smallestd) vNewToVar  =smalld*tend
          IF (vOldToVar  .LT.smallestd) vOldToVar  =smalld*tend

          d=((vNFrom-oned)*LOG(vOneBySq2Pi/SQRT(vNewFromVar))-(vNFrom-oned)/twod)- &
          &  (vNFrom*LOG(vOneBySq2Pi/SQRT(vOldFromVar))-(vNFrom)/twod)+            &
          &  ((vNTo+oned)*LOG(vOneBySq2Pi/SQRT(vNewToVar))-(vNTo+oned)/twod)-      &
          &  (vNTo*LOG(vOneBySq2Pi/SQRT(vOldToVar))-(vNTo)/twod)

          DTYPE(E_PCGaussian_EvaluateEnergyDifference)=-REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PCGaussian_EvaluateEnergyDifference)

        FUNCTION DTYPE(E_PCGaussian_EvaluateEnergyDifference_E_Merge)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,oldlabel_region,    &
        &        newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(E_PCGaussian)                                   :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(E_PCGaussian_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double), PARAMETER :: vOneBySq2Pi = oned/SQRT(twod*pid)
          REAL(ppm_kind_double)            :: vNFrom,vNTo
          REAL(ppm_kind_double)            :: vNewFromMean,vOldFromMean
          REAL(ppm_kind_double)            :: vNewToMean,vOldToMean
          REAL(ppm_kind_double)            :: vNewFromVar,vOldFromVar
          REAL(ppm_kind_double)            :: vNewToVar,vOldToVar
          REAL(ppm_kind_double)            :: d
          REAL(ppm_kind_double)            :: intensity
          REAL(MK)                         :: tmp

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          vNTo        =this%gCount(newlabel_region)
          vNFrom      =this%gCount(oldlabel_region)

          vNewFromMean=(this%gSums(oldlabel_region)-intensity)/(vNFrom - oned)
          vOldFromMean=this%gSums(oldlabel_region)/vNFrom
          vNewToMean  =(this%gSums(newlabel_region)+intensity)/(vNTo + oned)
          vOldToMean  =this%gSums(newlabel_region)/vNTo

          vNewFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region)- &
          &            intensity*intensity,vNewFromMean,vNFrom-oned)
          vOldFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region), &
          &            vOldFromMean,vNFrom)

          vNewToVar   =this%CalculateVariance(this%gSumsq(newlabel_region)+ &
          &            intensity*intensity,vNewToMean,vNTo+oned)
          vOldToVar   =this%CalculateVariance(this%gSumsq(newlabel_region), &
          &            vOldToMean,vNTo)

          ! it might happen that the variance is exactly 0 (if all pixel contain
          ! the same value). We therefore set it to a small value.
          IF (vNewFromVar.LT.smallestd) vNewFromVar=smalld*tend
          IF (vOldFromVar.LT.smallestd) vOldFromVar=smalld*tend
          IF (vNewToVar  .LT.smallestd) vNewToVar  =smalld*tend
          IF (vOldToVar  .LT.smallestd) vOldToVar  =smalld*tend

          tmp=this%CalculateKullbackLeiblerDistance(vNewToMean, &
          &   vOldFromMean,vNewToVar,vOldFromVar,vNTo+oned,vNFrom)
          e_merge=tmp.LT.this%RegionMergingThreshold

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PCGaussian_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ENDIF

          d=((vNFrom-oned)*LOG(vOneBySq2Pi/SQRT(vNewFromVar))-(vNFrom-oned)/twod)- &
          &  (vNFrom*LOG(vOneBySq2Pi/SQRT(vOldFromVar))-(vNFrom)/twod)+            &
          &  ((vNTo+oned)*LOG(vOneBySq2Pi/SQRT(vNewToVar))-(vNTo+oned)/twod)-      &
          &  (vNTo*LOG(vOneBySq2Pi/SQRT(vOldToVar))-(vNTo)/twod)

          DTYPE(E_PCGaussian_EvaluateEnergyDifference_E_Merge)=-REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PCGaussian_EvaluateEnergyDifference_E_Merge)

#if   __DIME == __3D
        SUBROUTINE E_PCPoisson_destroy(this,info)

          IMPLICIT NONE

          CLASS(E_PCPoisson)     :: this

          INTEGER, INTENT(  OUT) :: info

          info=0

        END SUBROUTINE E_PCPoisson_destroy
#endif

        FUNCTION DTYPE(E_PCPoisson_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,           &
        &        oldlabel_region,newlabel_region)

          IMPLICIT NONE

          CLASS(E_PCPoisson)                                    :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER      :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER      :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          REAL(MK)                                              :: DTYPE(E_PCPoisson_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vNFrom,vNTo
          REAL(ppm_kind_double) :: vNewFromMean,vOldFromMean
          REAL(ppm_kind_double) :: vNewToMean,vOldToMean
          REAL(ppm_kind_double) :: vNewFromVar,vOldFromVar
          REAL(ppm_kind_double) :: vNewToVar,vOldToVar
          REAL(ppm_kind_double) :: d
          REAL(ppm_kind_double) :: intensity

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PCPoisson_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          !this%gCount(oldlabel_region) can not be less than or equal to oned
          vNFrom      =this%gCount(oldlabel_region)
          !this%gCount(newlabel_region) can not be less than small
          vNTo        =this%gCount(newlabel_region)

          vNewFromMean=(this%gSums(oldlabel_region)-intensity)/(vNFrom - oned)
          vOldFromMean=this%gSums(oldlabel_region)/vNFrom

          vNewToMean  =(this%gSums(newlabel_region)+intensity)/(vNTo + oned)
          vOldToMean  =this%gSums(newlabel_region)/vNTo

          vNewFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region)- &
          &            intensity*intensity,vNewFromMean,vNFrom-oned)
          vOldFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region), &
          &            vOldFromMean,vNFrom)

          vNewToVar   =this%CalculateVariance(this%gSumsq(newlabel_region)+ &
          &            intensity*intensity,vNewToMean,vNTo+oned)
          vOldToVar   =this%CalculateVariance(this%gSumsq(newlabel_region), &
          &            vOldToMean,vNTo)

          ! it might happen that the variance is exactly 0 (if all pixel contain
          ! the same value). We therefore set it to a small value.
          IF (vNewFromVar.LT.smallestd) vNewFromVar=smalld*tend
          IF (vOldFromVar.LT.smallestd) vOldFromVar=smalld*tend
          IF (vNewToVar  .LT.smallestd) vNewToVar  =smalld*tend
          IF (vOldToVar  .LT.smallestd) vOldToVar  =smalld*tend

          d= LOG(vNewFromMean)*(this%gSums(oldlabel_region)-intensity)-(vNFrom-oned)*vNewFromMean- &
          & (LOG(vOldFromMean)*this%gSums(oldlabel_region)-vNFrom*vOldFromMean)+                   &
          &  LOG(vNewToMean)*(this%gSums(newlabel_region)+intensity)-(vNTo+oned)*vNewToMean-       &
          & (LOG(vOldToMean)*(this%gSums(newlabel_region))-vNTo*vOldToMean)

          DTYPE(E_PCPoisson_EvaluateEnergyDifference)=-REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PCPoisson_EvaluateEnergyDifference)

        FUNCTION DTYPE(E_PCPoisson_EvaluateEnergyDifference_E_Merge)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,oldlabel_region,   &
        &        newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(E_PCPoisson)                                    :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labels_
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labels_
#endif
          INTEGER,             DIMENSION(:),      INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: oldlabel
          INTEGER,                                INTENT(IN   ) :: newlabel
          INTEGER,                                INTENT(IN   ) :: oldlabel_region
          INTEGER,                                INTENT(IN   ) :: newlabel_region

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(E_PCPoisson_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vNFrom,vNTo
          REAL(ppm_kind_double) :: vNewFromMean,vOldFromMean
          REAL(ppm_kind_double) :: vNewToMean,vOldToMean
          REAL(ppm_kind_double) :: vNewFromVar,vOldFromVar
          REAL(ppm_kind_double) :: vNewToVar,vOldToVar
          REAL(ppm_kind_double) :: d
          REAL(ppm_kind_double) :: intensity
          REAL(MK)              :: tmp

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          vNFrom      =this%gCount(oldlabel_region)
          vNTo        =this%gCount(newlabel_region)

          vNewFromMean=(this%gSums(oldlabel_region)-intensity)/(vNFrom - oned)
          vOldFromMean=this%gSums(oldlabel_region)/vNFrom
          vNewToMean  =(this%gSums(newlabel_region)+intensity)/(vNTo + oned)
          vOldToMean  =this%gSums(newlabel_region)/vNTo

          vNewFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region)- &
          &            intensity*intensity,vNewFromMean,vNFrom-oned);
          vOldFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region), &
          &            vOldFromMean,vNFrom)

          vNewToVar   =this%CalculateVariance(this%gSumsq(newlabel_region)+ &
          &            intensity*intensity,vNewToMean,vNTo+oned)
          vOldToVar   =this%CalculateVariance(this%gSumsq(newlabel_region), &
          &            vOldToMean, vNTo)

          ! it might happen that the variance is exactly 0 (if all pixel contain
          ! the same value). We therefore set it to a small value.
          IF (vNewFromVar.LT.smallestd) vNewFromVar=smalld*tend
          IF (vOldFromVar.LT.smallestd) vOldFromVar=smalld*tend
          IF (vNewToVar  .LT.smallestd) vNewToVar  =smalld*tend
          IF (vOldToVar  .LT.smallestd) vOldToVar  =smalld*tend

          tmp=this%CalculateKullbackLeiblerDistance(vNewToMean, &
          &   vOldFromMean,vNewToVar,vOldFromVar,vNTo+oned,vNFrom)
          e_merge=tmp.LT.this%RegionMergingThreshold

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_PCPoisson_EvaluateEnergyDifference_E_Merge)=zero
             RETURN
          ENDIF

          d= LOG(vNewFromMean)*(this%gSums(oldlabel_region)-intensity)-(vNFrom-oned)*vNewFromMean- &
          & (LOG(vOldFromMean)*this%gSums(oldlabel_region)-vNFrom*vOldFromMean)+                   &
          &  LOG(vNewToMean)*(this%gSums(newlabel_region)+intensity)-(vNTo+oned)*vNewToMean-       &
          & (LOG(vOldToMean)*(this%gSums(newlabel_region))-vNTo*vOldToMean)

          DTYPE(E_PCPoisson_EvaluateEnergyDifference_E_Merge)=-REAL(this%m_Coefficient*d,MK)

        END FUNCTION DTYPE(E_PCPoisson_EvaluateEnergyDifference_E_Merge)
