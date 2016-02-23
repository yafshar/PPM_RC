#if   __DIME == __3D
        ! Constructor
        SUBROUTINE RCEnergyBaseClass_create(this,m_EnergyFunctional_, &
        &          m_Coefficient_,info,RegionMergingThreshold_)

          IMPLICIT NONE

          CLASS(RCEnergyBaseClass) :: this

          INTEGER,            INTENT(IN   ) :: m_EnergyFunctional_

          REAL(MK),           INTENT(IN   ) :: m_Coefficient_

          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: RegionMergingThreshold_

          this%m_Coefficient=REAL(m_Coefficient_,ppm_kind_double)
          this%m_EnergyFunctional=m_EnergyFunctional_
          IF (PRESENT(RegionMergingThreshold_)) THEN
             CALL this%Set(RegionMergingThreshold_)
          ENDIF

          info=0

        END SUBROUTINE RCEnergyBaseClass_create

        FUNCTION CalculateVariance(this,aSumSq,aMean,aN)
        !!! calculating an unbiased estimate of the population
        !!! variance from a finite sample of n observations
        !!!
        !!! Note
        !!! Because aSumSq and aN * aMean * aMean can be very similar numbers,
        !!! Cancellation can lead to the precision of the result to be much less
        !!! than the inherent precision of the floating-point arithmetic used to
        !!! perform the computation.
        !!! Thus this algorithm should not be used in practice.
        !!!
        !!! It can be changed to computing shifted data in future.
          IMPLICIT NONE

          CLASS(RCEnergyBaseClass)             :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: aSumSq
          REAL(ppm_kind_double), INTENT(IN   ) :: aMean
          REAL(ppm_kind_double), INTENT(IN   ) :: aN
          REAL(ppm_kind_double)                :: CalculateVariance

          IF (aN.LT.twod) THEN
             CalculateVariance=oned
          ELSE
             CalculateVariance=(aSumSq - aN * aMean * aMean)/(aN - oned)
          ENDIF

        END FUNCTION CalculateVariance

        FUNCTION CalculateScaledSphereVolume(this,RadiusX)
        !!! Get the sphere volume respecting the spacing of the label image.
        !!! A method commonly used by energies.
          IMPLICIT NONE

          CLASS(RCEnergyBaseClass) :: this

          REAL(MK),  INTENT(IN   ) :: RadiusX
          REAL(MK)                 :: CalculateScaledSphereVolume

          REAL(MK) :: RadiusY,RadiusZ

          RadiusY=RadiusX*pixel(1)/pixel(2)

          SELECT CASE (ppm_rc_dim)
          CASE (2)
             CalculateScaledSphereVolume = pi*RadiusX*RadiusY
          CASE (3)
             RadiusZ=RadiusX*pixel(1)/pixel(3)
             CalculateScaledSphereVolume = 4.0_MK/three*pi*RadiusX*RadiusY*RadiusZ
          END SELECT

        END FUNCTION CalculateScaledSphereVolume

        ! Constructor
        SUBROUTINE RCExternalEnergyBaseClass_Set(this,RegionMergingThreshold_)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          REAL(MK), INTENT(IN   ) :: RegionMergingThreshold_
          this%RegionMergingThreshold=RegionMergingThreshold_
        END SUBROUTINE RCExternalEnergyBaseClass_Set

        SUBROUTINE galloc_RCExternalEnergyBaseClass(this,Nregions_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: Nregions_
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          start_subroutine("galloc_RCExternalEnergyBaseClass")

          IF (Nregions_.LT.0) THEN
             fail("Nregions_ is less than 0!!!",ppm_error=ppm_error_fatal)
          ENDIF

          ALLOCATE(this%gCount(0:Nregions_), &
          &        this%gSums(0:Nregions_),  &
          &        this%gSumsq(0:Nregions_), &
          &        this%Rlabel(0:Nregions_), STAT=info)
          or_fail_alloc("gCount & gSums & gSumsq & Rlabel")

          !Initialize the variables
          this%Rlabel=-1
          this%gCount=zerod
          this%gSums =zerod
          this%gSumsq=zerod

          end_subroutine()

        END SUBROUTINE galloc_RCExternalEnergyBaseClass

        SUBROUTINE lalloc_RCExternalEnergyBaseClass(this,Nregions_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: Nregions_
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          start_subroutine("lalloc_RCExternalEnergyBaseClass")

          ALLOCATE(this%lCount(0:Nregions_), &
          &        this%lSumslSumsq(0:2*Nregions_+1),STAT=info)
          or_fail_alloc("lCount & lSumslSumsq")

          !Initialize the variables
          this%lCount     =zerod
          this%lSumslSumsq=zerod

          end_subroutine()

        END SUBROUTINE lalloc_RCExternalEnergyBaseClass

        SUBROUTINE grow_RCExternalEnergyBaseClass(this,Nregions_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: Nregions_
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          INTEGER :: nsize,i

          start_subroutine("grow_RCExternalEnergyBaseClass")

          nsize=SIZE(this%gCount)-1
          !upper bound of the array (0:nsize)

          IF (nsize.LT.Nregions_) THEN
             ALLOCATE(bufr(0:nsize),SOURCE=this%gCount(0:nsize),STAT=info)
             or_fail_alloc("bufr")
             DEALLOCATE(this%gCount,STAT=info)
             or_fail_dealloc("gCount")
             ALLOCATE(this%gCount(0:Nregions_),STAT=info)
             or_fail_alloc("gCount")

             FORALL (i=0:nsize)
                this%gCount(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:Nregions_)
                this%gCount(i)=zerod
             END FORALL

             bufr=this%gSums
             DEALLOCATE(this%gSums,STAT=info)
             or_fail_dealloc("gSums")
             ALLOCATE(this%gSums(0:Nregions_),STAT=info)
             or_fail_alloc("gSums")

             FORALL (i=0:nsize)
                this%gSums(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:Nregions_)
                this%gSums(i)=zerod
             END FORALL

             bufr=this%gSumsq
             DEALLOCATE(this%gSumsq,STAT=info)
             or_fail_dealloc("gSumsq")
             ALLOCATE(this%gSumsq(0:Nregions_),STAT=info)
             or_fail_alloc("gSumsq")

             FORALL (i=0:nsize)
                this%gSumsq(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:Nregions_)
                this%gSumsq(i)=zerod
             END FORALL

             DEALLOCATE(bufr,STAT=info)
             or_fail_dealloc("bufr")

             ALLOCATE(bufi(0:nsize),SOURCE=this%Rlabel(0:nsize),STAT=info)
             or_fail_alloc("bufi")
             DEALLOCATE(this%Rlabel,STAT=info)
             or_fail_dealloc("Rlabel")
             ALLOCATE(this%Rlabel(0:Nregions_),STAT=info)
             or_fail_alloc("Rlabel")

             FORALL (i=0:nsize)
                this%Rlabel(i)=bufi(i)
             END FORALL
             FORALL (i=nsize+1:Nregions_)
                this%Rlabel(i)=-1
             END FORALL

             DEALLOCATE(bufi,STAT=info)
             or_fail_dealloc("bufi")
          ENDIF

          end_subroutine()

        END SUBROUTINE grow_RCExternalEnergyBaseClass

        SUBROUTINE destroy_RCExternalEnergyBaseClass(this,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(  OUT) :: info

          start_subroutine("destroy_RCExternalEnergyBaseClass")

          IF (ALLOCATED(this%gCount)) THEN
             DEALLOCATE(this%gCount,STAT=info)
             or_fail_dealloc("gCount")
          ENDIF
          IF (ALLOCATED(this%gSums)) THEN
             DEALLOCATE(this%gSums,STAT=info)
             or_fail_dealloc("gSums")
          ENDIF
          IF (ALLOCATED(this%gSumsq)) THEN
             DEALLOCATE(this%gSumsq,STAT=info)
             or_fail_dealloc("gSumsq")
          ENDIF
          IF (ALLOCATED(this%Rlabel)) THEN
             DEALLOCATE(this%Rlabel,STAT=info)
             or_fail_dealloc("Rlabel")
          ENDIF
          IF (ALLOCATED(this%lCount)) THEN
             DEALLOCATE(this%lCount,STAT=info)
             or_fail_dealloc("lCount")
          ENDIF
          IF (ALLOCATED(this%lSumslSumsq)) THEN
             DEALLOCATE(this%lSumslSumsq,STAT=info)
             or_fail_dealloc("lSumslSumsq")
          ENDIF

          CALL this%destroy_(info)
          or_fail_dealloc("destroy_")

          end_subroutine()

        END SUBROUTINE destroy_RCExternalEnergyBaseClass
#endif

        FUNCTION DTYPE(RCExt_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,e_merge)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(RCExt_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: intensity

          INTEGER :: oldlabel_region,newlabel_region

          start_function("RCExt_EvaluateEnergyDifference")

          IF (newlabel.EQ.oldlabel) THEN
             DTYPE(RCExt_EvaluateEnergyDifference)=zero
             e_merge=.FALSE.
             GOTO 9999
          ENDIF

          ! translate region labels to region IDs using the hash table
          oldlabel_region=htable%search(oldlabel)
          newlabel_region=htable%search(newlabel)

          IF (oldlabel_region.EQ.htable_null.OR.newlabel_region.EQ.htable_null) THEN
             stdout("Energy can not be computed for a region who does not exist!")
             stdout(" Label1=",oldlabel," Label1 Index=",oldlabel_region)
             stdout(" Label2=",newlabel," Label2 Index=",newlabel_region)
             stdout(" at Coordinate =",coord)
             fail("The region does not exist!",ppm_error=ppm_error_fatal)
          ENDIF

          !The energy is defined to be -\inf in case the particle is the last
          !parent
          IF (this%gCount(oldlabel_region).LE.oneplusd) THEN
             DTYPE(RCExt_EvaluateEnergyDifference)=-bigs
             e_merge=.FALSE.
             GOTO 9999
          ENDIF

          IF (this%gCount(newlabel_region).LE.smalld) THEN
             DTYPE(RCExt_EvaluateEnergyDifference)=bigs
             e_merge=.FALSE.
             GOTO 9999
          ENDIF

          SELECT CASE (e_merge)
          CASE (.FALSE.)
             DTYPE(RCExt_EvaluateEnergyDifference) =     &
             & this%DTYPE(EvaluateEnergyDifference)_(    &
             & image_,labels_,coord,oldlabel,newlabel,   &
             & oldlabel_region,newlabel_region)

          CASE DEFAULT
             DTYPE(RCExt_EvaluateEnergyDifference) =     &
             & this%DTYPE(EvaluateEnergyDifference)__(   &
             & image_,labels_,coord,oldlabel,newlabel,   &
             & oldlabel_region,newlabel_region,e_merge)

          END SELECT

          IF (this%m_EnergyFunctional.GT.1000) THEN
             IF (oldlabel.EQ.0) THEN
#if   __DIME == __2D
                intensity=image_(coord(1),coord(2))
#elif __DIME == __3D
                intensity=image_(coord(1),coord(2),coord(3))
#endif
                IF (energy_coeff_balloon.GT.zero) THEN
                !!! outward flow
                   DTYPE(RCExt_EvaluateEnergyDifference)=   &
                   & DTYPE(RCExt_EvaluateEnergyDifference)- &
                   & energy_coeff_balloon*intensity
                ELSE
                   DTYPE(RCExt_EvaluateEnergyDifference)=   &
                   & DTYPE(RCExt_EvaluateEnergyDifference)+ &
                   & energy_coeff_balloon*(one-intensity)
                ENDIF
             ENDIF !oldlabel.EQ.0
          ENDIF !this%m_EnergyFunctional.GT.1000

          end_function()

        END FUNCTION DTYPE(RCExt_EvaluateEnergyDifference)


        FUNCTION DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)(this,   &
        &        image_,labels_,coord,oldlabel,newlabel,               &
        &        oldlabel_region,newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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

          REAL(MK)                                              :: DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vNewToMean,vOldFromMean,vNewToVar,vOldFromVar
          REAL(ppm_kind_double) :: vNTo,vNFrom
          REAL(ppm_kind_double) :: intensity

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          vNTo        =this%gCount(newlabel_region)+oned
          vNFrom      =this%gCount(oldlabel_region)

!           IF (vNTo+vNFrom.GT.10000._8) THEN
!              e_merge=.FALSE.
!              RETURN
!           ENDIF

          vNewToMean  =this%gSums(newlabel_region) + intensity
          vNewToMean  =vNewToMean/vNTo
          vOldFromMean=this%gSums(oldlabel_region)/vNFrom

          vNewToVar   =this%CalculateVariance(this%gSumsq(newlabel_region)+ &
          &            intensity*intensity,vNewToMean,vNTo)
          vOldFromVar =this%CalculateVariance(this%gSumsq(oldlabel_region), &
          &            vOldFromMean,vNFrom)

          DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)                  &
          &           =this%CalculateKullbackLeiblerDistance(vNewToMean, &
          &            vOldFromMean,vNewToVar,vOldFromVar,vNTo,vNFrom)

          e_merge=DTYPE(RCExt_EvaluateEnergyDifference_E_Merge).LT.this%RegionMergingThreshold

        END FUNCTION DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)

        FUNCTION DTYPE(CalculateTotalEnergy)(this,image_,labels_,Nm,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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

          REAL(MK)                                              :: DTYPE(CalculateTotalEnergy)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          DTYPE(CalculateTotalEnergy)=zero
          info=0
        END FUNCTION DTYPE(CalculateTotalEnergy)

#if   __DIME == __3D
        FUNCTION CalculateKullbackLeiblerDistance(this,aMu1,aMu2,aVar1_,aVar2_,aN1,aN2)
          !----------------------------------------------------------------------
          !  Kullback–Leibler divergence is a non-symmetric measure of the
          !  difference between two probability distributions.
          !----------------------------------------------------------------------

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)     :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: aMu1
          REAL(ppm_kind_double), INTENT(IN   ) :: aMu2
          REAL(ppm_kind_double), INTENT(IN   ) :: aVar1_
          REAL(ppm_kind_double), INTENT(IN   ) :: aVar2_
          REAL(ppm_kind_double), INTENT(IN   ) :: aN1
          REAL(ppm_kind_double), INTENT(IN   ) :: aN2
          REAL(MK)                             :: CalculateKullbackLeiblerDistance
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vMu12,vSumOfSq1,vSumOfSq2,vVar12
          REAL(ppm_kind_double) :: vDKL1,vDKL2,aVar1,aVar2

!           start_function("CalculateKullbackLeiblerDistance")

          vMu12=(aN1*aMu1+aN2*aMu2)/(aN1+aN2)

          !TOCHECK
          !I just added this because of accuracy
          aVar1=MERGE(smalld*tend,aVar1_,aVar1_.LT.smalld)
          aVar2=MERGE(smalld*tend,aVar2_,aVar2_.LT.smalld)

          vSumOfSq1=aVar1*(aN1-oned)+aN1*aMu1*aMu1
          vSumOfSq2=aVar2*(aN2-oned)+aN2*aMu2*aMu2

          vVar12=(oned/(aN1+aN2-oned))*(vSumOfSq1+vSumOfSq2-(aN1+aN2)*vMu12*vMu12)

          vDKL1=(aMu1-vMu12)*(aMu1-vMu12)/(twod*vVar12)+ &
          &     halfd*(aVar1/vVar12-oned-LOG(aVar1/vVar12))
          vDKL2=(aMu2-vMu12)*(aMu2-vMu12)/(twod*vVar12)+ &
                halfd*(aVar2/vVar12-oned-LOG(aVar2/vVar12))

          IF ((vDKL1+vDKL2).GT. bigsd.OR. &
          &   (vDKL1+vDKL2).LT.-bigsd) THEN
             CalculateKullbackLeiblerDistance=bigs
             RETURN
          ENDIF

          CalculateKullbackLeiblerDistance=REAL(vDKL1+vDKL2,MK)

!           end_function()

        END FUNCTION CalculateKullbackLeiblerDistance
#endif

        SUBROUTINE DTYPE(AddPoint)(this,image_,coord,label_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: label_
          INTEGER,                                INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: label_region,i

          start_subroutine("AddPoint")

          label_region=htable%search(label_)

          IF (label_region.EQ.htable_null) THEN
             fail("Adding point to the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          this%lCount(label_region)=this%lCount(label_region)+oned
          i=label_region*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  +intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)+intensity*intensity

          end_subroutine()

        END SUBROUTINE DTYPE(AddPoint)

#if   __DIME == __3D
        SUBROUTINE AddPoint_i(this,intensity_,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: intensity_
          INTEGER,           INTENT(IN   ) :: label_region_
          INTEGER,           INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: i

          start_subroutine("AddPoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Adding point to the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          intensity=REAL(intensity_,ppm_kind_double)

          this%lCount(label_region_)=this%lCount(label_region_)+oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  +intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)+intensity*intensity

          end_subroutine()

        END SUBROUTINE AddPoint_i

        SUBROUTINE AddPoint_rs(this,intensity_,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)     :: this

          REAL(ppm_kind_single), INTENT(IN   ) :: intensity_

          INTEGER,               INTENT(IN   ) :: label_region_
          INTEGER,               INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: i

          start_subroutine("AddPoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Adding point to the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          intensity=REAL(intensity_,ppm_kind_double)

          this%lCount(label_region_)=this%lCount(label_region_)+oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  +intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)+intensity*intensity

          end_subroutine()

        END SUBROUTINE AddPoint_rs

        SUBROUTINE AddPoint_r(this,intensity,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)     :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: intensity

          INTEGER,               INTENT(IN   ) :: label_region_
          INTEGER,               INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: i

          start_subroutine("AddPoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Adding point to the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          this%lCount(label_region_)=this%lCount(label_region_)+oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  +intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)+intensity*intensity

          end_subroutine()

        END SUBROUTINE AddPoint_r
#endif

        SUBROUTINE DTYPE(RemovePoint)(this,image_,coord,label_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: label_
          INTEGER,                                INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: label_region,i

          start_subroutine("RemovePoint")

          label_region=htable%search(label_)

          IF (label_region.EQ.htable_null) THEN
             fail("Removing point from the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          this%lCount(label_region)=this%lCount(label_region)-oned
          i=label_region*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  -intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)-intensity*intensity

          end_subroutine()

        END SUBROUTINE DTYPE(RemovePoint)

#if   __DIME == __3D
        SUBROUTINE RemovePoint_i(this,intensity_,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: intensity_
          INTEGER,           INTENT(IN   ) :: label_region_
          INTEGER,           INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: i

          start_subroutine("RemovePoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Removing point from the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          intensity=REAL(intensity_,ppm_kind_double)

          this%lCount(label_region_)=this%lCount(label_region_)-oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  -intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)-intensity*intensity

          end_subroutine()

        END SUBROUTINE RemovePoint_i

        SUBROUTINE RemovePoint_rs(this,intensity_,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)     :: this

          REAL(ppm_kind_single), INTENT(IN   ) :: intensity_

          INTEGER,               INTENT(IN   ) :: label_region_
          INTEGER,               INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: i

          start_subroutine("RemovePoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Removing point from the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          intensity=REAL(intensity_,ppm_kind_double)

          this%lCount(label_region_)=this%lCount(label_region_)-oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  -intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)-intensity*intensity

          end_subroutine()

        END SUBROUTINE RemovePoint_rs

        SUBROUTINE RemovePoint_r(this,intensity,label_region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)     :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: intensity

          INTEGER,               INTENT(IN   ) :: label_region_
          INTEGER,               INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: i

          start_subroutine("RemovePoint")

          IF (label_region_.EQ.htable_null) THEN
             fail("Removing point from the region which does not exist", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          this%lCount(label_region_)=this%lCount(label_region_)-oned
          i=label_region_*2
          this%lSumslSumsq(i)  =this%lSumslSumsq(i)  -intensity
          this%lSumslSumsq(i+1)=this%lSumslSumsq(i+1)-intensity*intensity

          end_subroutine()

        END SUBROUTINE RemovePoint_r

        SUBROUTINE KillRegion(this,region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: region_
          INTEGER,           INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64) :: key

          start_subroutine("KillRegion")

          SELECT CASE (region_)
          CASE (0)
             this%gCount(0)=zerod
             this%gSums(0) =zerod
             this%gSumsq(0)=zerod

          CASE DEFAULT
             this%gCount(region_)=smalld
             key=INT(this%Rlabel(region_),ppm_kind_int64)
             IF (htable%search(key).NE.htable_null) THEN
                CALL htable%remove(key,info,.TRUE.)
                or_fail("hash remove")
             ENDIF
             this%Rlabel(region_)=-1
             this%gSums(region_) =zerod
             this%gSumsq(region_)=zerod

          END SELECT

          end_subroutine()

        END SUBROUTINE KillRegion

        ! Constructor
        SUBROUTINE RCInternalEnergyBaseClass_Set(this,RegionMergingThreshold_)

          IMPLICIT NONE

          CLASS(RCInternalEnergyBaseClass) :: this

          REAL(MK), INTENT(IN   ) :: RegionMergingThreshold_
        END SUBROUTINE RCInternalEnergyBaseClass_Set
#endif

        FUNCTION DTYPE(RCInt_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel,e_merge)

          IMPLICIT NONE

          CLASS(RCInternalEnergyBaseClass)                     :: this

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

          LOGICAL,                                INTENT(INOUT) :: e_merge

          REAL(MK)                                              :: DTYPE(RCInt_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          IF (newlabel.EQ.oldlabel) THEN
             DTYPE(RCInt_EvaluateEnergyDifference) = zero
             GOTO 9999
          ENDIF

          DTYPE(RCInt_EvaluateEnergyDifference) =         &
          & this%DTYPE(EvaluateEnergyDifference)_(image_, &
          & labels_,coord,oldlabel,newlabel)

          IF (this%m_EnergyFunctional.GT.1000) THEN
             IF (oldlabel.EQ.0) THEN
                DTYPE(RCInt_EvaluateEnergyDifference)   = &
                & DTYPE(RCInt_EvaluateEnergyDifference) - &
                & energy_coeff_outward_flow
             ELSE IF (newlabel.EQ.0) THEN
                DTYPE(RCInt_EvaluateEnergyDifference)   = &
                & DTYPE(RCInt_EvaluateEnergyDifference) + &
                & energy_coeff_outward_flow
             ENDIF
          ENDIF

        9999 CONTINUE
          RETURN
        END FUNCTION DTYPE(RCInt_EvaluateEnergyDifference)

        SUBROUTINE DTYPE(UpdateStatisticsWhenJump)(this, &
        &          image_,coord,aFromLabel,aToLabel,info)
          !Update the statistics of the propagating and the loser region.
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
          CLASS(RCExternalEnergyBaseClass)                      :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: image_
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
#endif

          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,                                INTENT(IN   ) :: aFromLabel
          INTEGER,                                INTENT(IN   ) :: aToLabel
          INTEGER,                                INTENT(  OUT) :: info
          !---------------------------------------------------------------------
          !  Local variables
          !---------------------------------------------------------------------
          REAL(ppm_kind_double) :: intensity

          INTEGER :: label_region

          start_subroutine("UpdateStatisticsWhenJump")
          ! The default behaviour is to call AddPoint(...) with the
          ! label where the point will move to and RemovePoint(...)
          ! with the label of which the point is "coming" from.
          ! For some energies it might be beneficial to overwrite
          ! this method.
#if   __DIME == __2D
          intensity=REAL(image_(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(image_(coord(1),coord(2),coord(3)),ppm_kind_double)
#endif

          label_region=htable%search(aToLabel)

          check_true(<#label_region.NE.htable_null#>, &
          & "Fail!!!, There should be an aToLabel available.")

          CALL this%AddPoint(intensity,label_region,info)
          or_fail("this%AddPoint")

          label_region=htable%search(aFromLabel)

          check_true(<#label_region.NE.htable_null#>, &
          & "Fail!!!, There should be an aFromLabel available.")

          CALL this%RemovePoint(intensity,label_region,info)
          or_fail("this%RemovePoint")

          end_subroutine()

        END SUBROUTINE DTYPE(UpdateStatisticsWhenJump)

#if   __DIME == __3D
        SUBROUTINE UpdateStatistics(this,info,sendrecv)
          !Update the statistics of the propagating and the loser region.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(  OUT) :: info

          LOGICAL, OPTIONAL, INTENT(IN   ) :: sendrecv
          !---------------------------------------------------------------------
          !  Local variables
          !---------------------------------------------------------------------
          INTEGER :: nsize,nsize_
          INTEGER :: i

          LOGICAL :: sendrecv1,sendrecv2

          start_subroutine("UpdateStatistics")

          IF (PRESENT(sendrecv)) THEN
             sendrecv1=sendrecv
             sendrecv2=.NOT.sendrecv
          ELSE
             sendrecv1=.TRUE.
             sendrecv2=.TRUE.
          ENDIF

          nsize=SIZE(this%lCount)

#ifdef __MPI
          IF (sendrecv1) THEN
             CALL MPI_Iallreduce(MPI_IN_PLACE,this%lCount(0), &
             &    nsize,MPI_DOUBLE_PRECISION,MPI_SUM,comm,requestCount,info)
             or_fail_MPI("MPI_Iallreduce")

             CALL MPI_Iallreduce(MPI_IN_PLACE,this%lSumslSumsq(0), &
             &    nsize*2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,requestSums,info)
             or_fail_MPI("MPI_Iallreduce")
          ENDIF !sendrecv1
#endif

          IF (sendrecv2) THEN
             nsize_=nsize-1

#ifdef __MPI
             CALL MPI_Wait(requestCount,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
#endif

             FORALL (i=0:nsize_,ABS(this%lCount(i)).GT.halfd)
                this%gCount(i)=this%gCount(i)+this%lCount(i)
             END FORALL

             DEALLOCATE(this%lCount,STAT=info)
             or_fail_dealloc("local array deallocation failed!")

             DO i=0,nsize_
                IF (this%gCount(i).LE.smallestd) THEN
                   CALL this%KillRegion(i,info)
                   or_fail("KillRegion")
                ENDIF
             ENDDO

#ifdef __MPI
             CALL MPI_Wait(requestSums,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
#endif

             FORALL (i=0:nsize_,this%gCount(i).GE.oned)
                this%gSums(i) =this%gSums(i) +this%lSumslSumsq(2*i)
                this%gSumsq(i)=this%gSumsq(i)+this%lSumslSumsq(2*i+1)
             END FORALL

             DEALLOCATE(this%lSumslSumsq,STAT=info)
             or_fail_dealloc("local array deallocation failed!")
          ENDIF !sendrecv2

          end_subroutine()

        END SUBROUTINE UpdateStatistics

        SUBROUTINE RemoveNotSignificantRegions(this,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER ::i

          start_subroutine("RemoveNotSignificantRegions")

           DO i=1,SIZE(e_data%gCount)-1
              IF (e_data%gCount(i).LE.RegionSizeThreshold.AND.e_data%gCount(i).GT.oneminusd) THEN
                 CALL e_data%RemoveFGRegion(i,info)
                 or_fail("RemoveFGRegion")
              ENDIF
           ENDDO

          end_subroutine()

        END SUBROUTINE RemoveNotSignificantRegions

        SUBROUTINE RemoveFGRegion(this,region_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: region_
          INTEGER,           INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpl_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpl_3d
          INTEGER,             DIMENSION(:),     POINTER :: Nm
          INTEGER                                        :: label_region
          INTEGER                                        :: i,j,k

          start_subroutine("RemoveFGRegion")

          label_region=this%Rlabel(region_)

          NULLIFY(wpl_2d,wpl_3d)

          sbpitr => mesh%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             SELECT CASE (SIZE(Nm))
             CASE (2)
                CALL sbpitr%get_field(labels,wpl_2d,info)
                or_fail("Failed to get field wpl_2d data.")

                FORALL (i=1:Nm(1),j=1:Nm(2),ABS(wpl_2d(i,j)).EQ.label_region)
                   wpl_2d(i,j)=0
                END FORALL

             CASE (3)
                CALL sbpitr%get_field(labels,wpl_3d,info)
                or_fail("Failed to get field wpl_2d data.")

                FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3),ABS(wpl_3d(i,j,k)).EQ.label_region)
                   wpl_3d(i,j,k)=0
                END FORALL

             END SELECT
             sbpitr => mesh%subpatch%next()
          ENDDO !WHILE (ASSOCIATED(sbpitr))

          this%gCount(0)=this%gCount(0)+this%gCount(region_)
          this%gSums(0) =this%gSums(0) +this%gSums(region_)
          this%gSumsq(0)=this%gSumsq(0)+this%gSumsq(region_)

          this%gCount(region_)=zerod

          end_subroutine()

        END SUBROUTINE RemoveFGRegion
#endif

        FUNCTION DTYPE(E_Gamma_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel)

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          CLASS(E_Gamma)                                        :: this

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

          REAL(MK)                                              :: DTYPE(E_Gamma_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER                            :: nchgold,nchgnew,i
          INTEGER, DIMENSION(__DIME)         :: ld
#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
#endif

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_Gamma_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

          !!! the bounds of tmplabels(1:3,1:3,1:3)
          tmplabels => &
#if   __DIME == __2D
          & labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)
#elif __DIME == __3D
          & labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)
#endif

          nchgold = 0
          nchgnew = 0
          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ld=FG_ConnectivityType%NeighborsPoints(:,i)+2

#if   __DIME == __2D
             IF (ABS(tmplabels(ld(1),ld(2))).NE.oldlabel) nchgold=nchgold+1
             IF (ABS(tmplabels(ld(1),ld(2))).NE.newlabel) nchgnew=nchgnew+1
#elif __DIME == __3D
             IF (ABS(tmplabels(ld(1),ld(2),ld(3))).NE.oldlabel) nchgold=nchgold+1
             IF (ABS(tmplabels(ld(1),ld(2),ld(3))).NE.newlabel) nchgnew=nchgnew+1
#endif
          ENDDO

          DTYPE(E_Gamma_EvaluateEnergyDifference)= &
          & REAL(this%m_Coefficient,MK)*REAL(nchgnew-nchgold,MK)

        END FUNCTION DTYPE(E_Gamma_EvaluateEnergyDifference)

#if   __DIME == __3D
        SUBROUTINE E_ContourLengthApprox_destroy(this,info)
        !!! This subroutine will destroy the allocated MASK array

          IMPLICIT NONE

          CLASS(E_ContourLengthApprox) :: this

          INTEGER,       INTENT(  OUT) :: info

          start_subroutine("E_ContourLengthApprox_destroy")

          DEALLOCATE(this%NeighborsPoints,STAT=info)
          or_fail_dealloc("this%NeighborsPoints")

          end_subroutine()

        END SUBROUTINE E_ContourLengthApprox_destroy
#endif

        SUBROUTINE DTYPE(E_ContourLengthApprox_PrepareEnergyCalculation)(this,info)

          IMPLICIT NONE

          CLASS(E_ContourLengthApprox) :: this

          INTEGER,       INTENT(  OUT) :: info

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
          !!! get the curvature radius in pixel
          e_lX=CEILING(this%m_Radius)
          !!! scaling factor
          sX=REAL(e_lX*e_lX,MK)

          e_lY=CEILING(this%m_Radius*pixel(1)/pixel(2))
          sY=REAL(e_lY*e_lY,MK)

          !!! this index is used for index shift in the MASK array
          lX=e_lX+1
          lY=e_lY+1

          !!! the size of the MASK array
          sizeX=e_lX*2+1
          sizeY=e_lY*2+1

          !!! the Area of the circle in 2D or the sphere in 3D
          this%m_Volume=this%CalculateScaledSphereVolume(this%m_Radius)

          l=0

#if   __DIME == __2D
          this%m_Prefactor=three*pi/this%m_Radius

          DO j=-e_lY,e_lY
             DO i=-e_lX,e_lX
                s=REAL(i*i)/sX+REAL(j*j)/sY
                IF (s.LE.oneplus) THEN
                   l=l+1
                ENDIF
             ENDDO
          ENDDO
          !!! Filling the MASK array based on the hyper sphere location

          this%NumberOfNeighbors=l

          ALLOCATE(this%NeighborsPoints(__DIME,l),STAT=info)
          or_fail_alloc("this%NeighborsPoints")

          l=0
          DO j=-e_lY,e_lY
             DO i=-e_lX,e_lX
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
          this%m_Prefactor=16.0_MK/(three*this%m_Radius)
          !!! The constant that depends on the dimension d and on Radius

          e_lZ=CEILING(this%m_Radius*pixel(1)/pixel(3))
          sZ=REAL(e_lZ*e_lZ,MK)

          !!! this index is used for index shift in the MASK array
          lZ=e_lZ+1

          !!! the size of the MASK array
          sizeZ=e_lZ*2+1

          DO k=-e_lZ,e_lZ
             DO j=-e_lY,e_lY
                DO i=-e_lX,e_lX
                   s=REAL(i*i)/sX+(j*j)/sY+(k*k)/sZ
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
          DO k=-e_lZ,e_lZ
             DO j=-e_lY,e_lY
                DO i=-e_lX,e_lX
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
          !!! Filling the NeighborsPoints array based on the hyper sphere location,
#endif
          end_subroutine()

        END SUBROUTINE DTYPE(E_ContourLengthApprox_PrepareEnergyCalculation)


        FUNCTION DTYPE(E_ContourLengthApprox_EvaluateEnergyDifference)(this, &
        &        image_,labels_,coord,oldlabel,newlabel)

          IMPLICIT NONE

          CLASS(E_ContourLengthApprox)                          :: this

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

          REAL(MK)                                              :: DTYPE(E_ContourLengthApprox_EvaluateEnergyDifference)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: vCurvatureFlow
          REAL(MK) :: vNTo,vNFrom

          INTEGER                            :: i,j,l
#if   __DIME == __2D
          INTEGER, DIMENSION(:,:),   POINTER :: tmplabels
#elif __DIME == __3D
          INTEGER                            :: k
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
#endif

          IF (this%m_Coefficient.LT.smallestd) THEN
             DTYPE(E_ContourLengthApprox_EvaluateEnergyDifference)=zero
             RETURN
          ENDIF

#if   __DIME == __2D
          tmplabels => labels_(coord(1)-e_lX:coord(1)+e_lX, &
          &                    coord(2)-e_lY:coord(2)+e_lY)
#elif __DIME == __3D
          tmplabels => labels_(coord(1)-e_lX:coord(1)+e_lX, &
          &                    coord(2)-e_lY:coord(2)+e_lY, &
          &                    coord(3)-e_lZ:coord(3)+e_lZ)
#endif

          IF (oldlabel.EQ.0) THEN !growing
             ! This is a point on the contour (innerlist) OR
             ! touching the contour (Outer list)
             vNTo=zero

             DO l=1,this%NumberOfNeighbors
                i=this%NeighborsPoints(1,l)
                j=this%NeighborsPoints(2,l)
#if   __DIME == __2D
                IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
#elif __DIME == __3D
                k=this%NeighborsPoints(3,l)
                IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
#endif
                   vNTo=vNTo+one
                ENDIF
             ENDDO

             vCurvatureFlow=-(vNTo/this%m_Volume-half)
          ELSE IF (newlabel.EQ.0) THEN !proper shrinking
             ! This is a point on the contour (innerlist) OR
             ! touching the contour (Outer list)
             vNFrom=zero

             DO l=1,this%NumberOfNeighbors
                i=this%NeighborsPoints(1,l)
                j=this%NeighborsPoints(2,l)
#if   __DIME == __2D
                IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
#elif __DIME == __3D
                k=this%NeighborsPoints(3,l)
                IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
#endif
                   vNFrom=vNFrom+one
                ENDIF
             ENDDO

             vCurvatureFlow=vNFrom/this%m_Volume-half
          ELSE ! fighting fronts
             vNFrom=zero
             vNTo  =zero

             DO l=1,this%NumberOfNeighbors
                i=this%NeighborsPoints(1,l)
                j=this%NeighborsPoints(2,l)
#if   __DIME == __2D
                IF (ABS(tmplabels(i,j)).EQ.oldlabel) THEN
                   vNFrom=vNFrom+one
                ELSE IF (ABS(tmplabels(i,j)).EQ.newlabel) THEN
#elif __DIME == __3D
                k=this%NeighborsPoints(3,l)
                IF (ABS(tmplabels(i,j,k)).EQ.oldlabel) THEN
                   vNFrom=vNFrom+one
                ELSE IF (ABS(tmplabels(i,j,k)).EQ.newlabel) THEN
#endif
                   vNTo=vNTo+one
                ENDIF
             ENDDO

!              vCurvatureFlow=-(vNTo/this%m_Volume-half)+vNFrom/this%m_Volume-half
             vCurvatureFlow=(-vNTo+vNFrom)/this%m_Volume
          ENDIF !oldlabel.EQ.0

          DTYPE(E_ContourLengthApprox_EvaluateEnergyDifference)= &
          & REAL(this%m_Coefficient,MK)*this%m_Prefactor*vCurvatureFlow

        END FUNCTION DTYPE(E_ContourLengthApprox_EvaluateEnergyDifference)

