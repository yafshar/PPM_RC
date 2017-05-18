#if   __DIME == __3D
        ! Constructor
        SUBROUTINE RCEnergyBaseClass_create(this,m_EnergyFunctionalIn, &
        &          m_Coefficient_,info,RegionMergingThresholdIn)

          IMPLICIT NONE

          CLASS(RCEnergyBaseClass) :: this

          INTEGER,            INTENT(IN   ) :: m_EnergyFunctionalIn

          REAL(MK),           INTENT(IN   ) :: m_Coefficient_

          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: RegionMergingThresholdIn

          this%m_Coefficient=REAL(m_Coefficient_,ppm_kind_double)
          this%m_EnergyFunctional=m_EnergyFunctionalIn
          IF (PRESENT(RegionMergingThresholdIn)) THEN
             CALL this%Set(RegionMergingThresholdIn)
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
        SUBROUTINE RCExternalEnergyBaseClass_Set(this,RegionMergingThresholdIn)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          REAL(MK), INTENT(IN   ) :: RegionMergingThresholdIn
          this%RegionMergingThreshold=RegionMergingThresholdIn
        END SUBROUTINE RCExternalEnergyBaseClass_Set

        SUBROUTINE galloc_RCExternalEnergyBaseClass(this,NregionsIn,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: NregionsIn
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          start_subroutine("galloc_RCExternalEnergyBaseClass")

          IF (NregionsIn.LT.0) THEN
             fail("NregionsIn is less than 0!!!",ppm_error=ppm_error_fatal)
          ENDIF

          !----------------------------------------------------------------------
          ! Allocate and Initialize the variables
          !----------------------------------------------------------------------
          ALLOCATE(this%gCount(0:NregionsIn),this%gSums(0:NregionsIn), &
          &        this%gSumsq(0:NregionsIn),SOURCE=zerod,STAT=info)
          or_fail_alloc("gCount, gSums & gSumsq")

          ALLOCATE(this%Rlabel(0:NregionsIn),SOURCE=-1,STAT=info)
          or_fail_alloc("Rlabel")

          end_subroutine()

        END SUBROUTINE galloc_RCExternalEnergyBaseClass

        SUBROUTINE lalloc_RCExternalEnergyBaseClass(this,NregionsIn,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: NregionsIn
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          start_subroutine("lalloc_RCExternalEnergyBaseClass")

          !----------------------------------------------------------------------
          ! Allocate and Initialize the variables
          !----------------------------------------------------------------------
          ALLOCATE(this%lCount(0:NregionsIn),this%lSumslSumsq(0:2*NregionsIn+1),SOURCE=zerod,STAT=info)
          or_fail_alloc("lCount & lSumslSumsq")

          end_subroutine()

        END SUBROUTINE lalloc_RCExternalEnergyBaseClass

        SUBROUTINE grow_RCExternalEnergyBaseClass(this,NregionsIn,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(IN   ) :: NregionsIn
          !!! number of regions except background
          INTEGER,           INTENT(  OUT) :: info

          INTEGER :: nsize,i

          start_subroutine("grow_RCExternalEnergyBaseClass")

          nsize=SIZE(this%gCount)-1
          !upper bound of the array (0:nsize)

          IF (nsize.LT.NregionsIn) THEN
             ALLOCATE(bufr(0:nsize),SOURCE=this%gCount(0:nsize),STAT=info)
             or_fail_alloc("bufr")
             DEALLOCATE(this%gCount,STAT=info)
             or_fail_dealloc("gCount")
             ALLOCATE(this%gCount(0:NregionsIn),STAT=info)
             or_fail_alloc("gCount")

             FORALL (i=0:nsize)
                this%gCount(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:NregionsIn)
                this%gCount(i)=zerod
             END FORALL

             bufr=this%gSums
             DEALLOCATE(this%gSums,STAT=info)
             or_fail_dealloc("gSums")
             ALLOCATE(this%gSums(0:NregionsIn),STAT=info)
             or_fail_alloc("gSums")

             FORALL (i=0:nsize)
                this%gSums(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:NregionsIn)
                this%gSums(i)=zerod
             END FORALL

             bufr=this%gSumsq
             DEALLOCATE(this%gSumsq,STAT=info)
             or_fail_dealloc("gSumsq")
             ALLOCATE(this%gSumsq(0:NregionsIn),STAT=info)
             or_fail_alloc("gSumsq")

             FORALL (i=0:nsize)
                this%gSumsq(i)=bufr(i)
             END FORALL
             FORALL (i=nsize+1:NregionsIn)
                this%gSumsq(i)=zerod
             END FORALL

             DEALLOCATE(bufr,STAT=info)
             or_fail_dealloc("bufr")

             ALLOCATE(bufi(0:nsize),SOURCE=this%Rlabel(0:nsize),STAT=info)
             or_fail_alloc("bufi")
             DEALLOCATE(this%Rlabel,STAT=info)
             or_fail_dealloc("Rlabel")
             ALLOCATE(this%Rlabel(0:NregionsIn),STAT=info)
             or_fail_alloc("Rlabel")

             FORALL (i=0:nsize)
                this%Rlabel(i)=bufi(i)
             END FORALL
             FORALL (i=nsize+1:NregionsIn)
                this%Rlabel(i)=-1
             END FORALL

             DEALLOCATE(bufi,STAT=info)
             or_fail_dealloc("bufi")
          ENDIF

          end_subroutine()

        END SUBROUTINE grow_RCExternalEnergyBaseClass

        SUBROUTINE shrink_RCExternalEnergyBaseClass(this,info,shrinkage_ratio)
          !shrink the array and sort them according to the region labels
          USE ppm_module_util_qsort, ONLY : ppm_util_qsort
          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass) :: this

          INTEGER,           INTENT(  OUT) :: info
          INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
          !!! OPTIONAL shrinkage_ratio (positive value).
          !!! If the size of hash table is shrinkage_ratio times bigger than the
          !!! real elements inside table, we reduce the table size

          INTEGER                 :: NregionsT
          INTEGER                 :: nsize,i,j
          INTEGER                 :: shrinkage_ratio_
          INTEGER(ppm_kind_int64) :: li

          start_subroutine("shrink_RCExternalEnergyBaseClass")

          shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
          shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

          NregionsT=this%size(bufl)

          nsize=SIZE(bufl)
          !upper bound of the array (0:nsize)

          IF (nsize.GT.shrinkage_ratio_*NregionsT) THEN
             CALL MOVE_ALLOC(this%gCount,bufr)
             ALLOCATE(this%gCount(0:NregionsT),STAT=info)
             or_fail_alloc("gCount")

             this%gCount(0)=bufr(0)
             j=0
             DO i=1,nsize
                IF (bufl(i)) THEN
                   j=j+1
                   this%gCount(j)=bufr(i)
                ENDIF
             ENDDO

             CALL MOVE_ALLOC(this%gSums,bufr)
             ALLOCATE(this%gSums(0:NregionsT),STAT=info)
             or_fail_alloc("gSums")

             this%gSums(0)=bufr(0)
             j=0
             DO i=1,nsize
                IF (bufl(i)) THEN
                   j=j+1
                   this%gSums(j)=bufr(i)
                ENDIF
             ENDDO

             CALL MOVE_ALLOC(this%gSumsq,bufr)
             ALLOCATE(this%gSumsq(0:NregionsT),STAT=info)
             or_fail_alloc("gSumsq")

             this%gSumsq(0)=bufr(0)
             j=0
             DO i=1,nsize
                IF (bufl(i)) THEN
                   j=j+1
                   this%gSumsq(j)=bufr(i)
                ENDIF
             ENDDO

             DEALLOCATE(bufr,STAT=info)
             or_fail_dealloc("bufr")

             CALL MOVE_ALLOC(this%Rlabel,bufi)
             ALLOCATE(this%Rlabel(0:NregionsT),STAT=info)
             or_fail_alloc("Rlabel")

             this%Rlabel(0)=bufi(0)
             j=0
             DO i=1,nsize
                IF (bufl(i)) THEN
                   j=j+1
                   this%Rlabel(j)=bufi(i)
                ENDIF
             ENDDO

             DEALLOCATE(bufi,bufl,STAT=info)
             or_fail_dealloc("bufi & bufl")

             NULLIFY(bufs)
             CALL ppm_util_qsort(this%Rlabel,bufs,info)
             or_fail("ppm_util_qsort")

             ! this%Rlabel starts from dimension 0 which inside the sorting routine
             ! will be as 1, so I need to reduce one for correct indexing
             bufs=bufs-1

             ALLOCATE(bufr(0:NregionsT),STAT=info)
             or_fail_alloc("bufr")

             DO i=1,NregionsT+1
                bufr(i-1)=this%gCount(bufs(i))
             ENDDO
             FORALL (i=0:NregionsT)
                this%gCount(i)=bufr(i)
             END FORALL

             DO i=1,NregionsT+1
                bufr(i-1)=this%gSums(bufs(i))
             ENDDO
             FORALL (i=0:NregionsT)
                this%gSums(i)=bufr(i)
             END FORALL

             DO i=1,NregionsT+1
                bufr(i-1)=this%gSumsq(bufs(i))
             ENDDO
             FORALL (i=0:NregionsT)
                this%gSumsq(i)=bufr(i)
             END FORALL

             DEALLOCATE(bufr,STAT=info)
             or_fail_dealloc("bufr")

             ALLOCATE(bufi(0:NregionsT),STAT=info)
             or_fail_alloc("bufi")

             DO i=1,NregionsT+1
                bufi(i-1)=this%Rlabel(bufs(i))
             ENDDO
             FORALL (i=0:NregionsT)
                this%Rlabel(i)=bufi(i)
             END FORALL

             DEALLOCATE(bufi,bufs,STAT=info)
             or_fail_dealloc("bufi & bufs")
             NULLIFY(bufs)

             CALL htable%destroy(info)
             or_fail("Failed to destroy the hash table!")

             CALL htable%create(NregionsT,info)
             or_fail("create htable")

             DO i=0,NregionsT
                li=INT(this%Rlabel(i),ppm_kind_int64)
                CALL htable%insert(li,i,info)
                or_fail("hash insert")
             ENDDO
          ELSE
             DEALLOCATE(bufl,STAT=info)
             or_fail_dealloc("bufl")
          ENDIF

          end_subroutine()

        END SUBROUTINE shrink_RCExternalEnergyBaseClass

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

        INTEGER FUNCTION size_RCExternalEnergyBaseClass(this,MASK)
        !!! Returns the number of connected FG regions without Background
          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                            :: this

          LOGICAL, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: MASK

          INTEGER :: nsize,i

          IF (PRESENT(MASK)) THEN
             nsize=SIZE(e_data%gCount)-1
             IF (ALLOCATED(MASK)) DEALLOCATE(MASK)
             ALLOCATE(MASK(nsize),SOURCE=.FALSE.)
             FORALL (i=1:nsize,e_data%gCount(i).GT.halfd)
                MASK(i)=.TRUE.
             END FORALL
             size_RCExternalEnergyBaseClass=COUNT(MASK)
          ELSE
             size_RCExternalEnergyBaseClass=COUNT(e_data%gCount.GT.halfd)-1
          ENDIF
        END FUNCTION size_RCExternalEnergyBaseClass
#endif

        FUNCTION DTYPE(RCExt_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,e_merge)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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
             & imageIn,labelsIn,coord,oldlabel,newlabel, &
             & oldlabel_region,newlabel_region)

          CASE DEFAULT
             DTYPE(RCExt_EvaluateEnergyDifference) =     &
             & this%DTYPE(EvaluateEnergyDifference)__(   &
             & imageIn,labelsIn,coord,oldlabel,newlabel, &
             & oldlabel_region,newlabel_region,e_merge)

          END SELECT

          IF (this%m_EnergyFunctional.GT.1000) THEN
             IF (oldlabel.EQ.0) THEN
#if   __DIME == __2D
                intensity=imageIn(coord(1),coord(2))
#elif __DIME == __3D
                intensity=imageIn(coord(1),coord(2),coord(3))
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


        FUNCTION DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,           &
        &        oldlabel_region,newlabel_region,e_merge)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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

          REAL(MK)                                              :: DTYPE(RCExt_EvaluateEnergyDifference_E_Merge)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: vNewToMean,vOldFromMean,vNewToVar,vOldFromVar
          REAL(ppm_kind_double) :: vNTo,vNFrom
          REAL(ppm_kind_double) :: intensity

#if   __DIME == __2D
          intensity=REAL(imageIn(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(imageIn(coord(1),coord(2),coord(3)),ppm_kind_double)
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

        FUNCTION DTYPE(CalculateTotalEnergy)(this,imageIn,labelsIn,Nm,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

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

        SUBROUTINE DTYPE(AddPoint)(this,imageIn,coord,label_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
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
          intensity=REAL(imageIn(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(imageIn(coord(1),coord(2),coord(3)),ppm_kind_double)
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

        SUBROUTINE DTYPE(RemovePoint)(this,imageIn,coord,label_,info)

          IMPLICIT NONE

          CLASS(RCExternalEnergyBaseClass)                      :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
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
          intensity=REAL(imageIn(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(imageIn(coord(1),coord(2),coord(3)),ppm_kind_double)
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
          !!! Region ID (not region label!)
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
        SUBROUTINE RCInternalEnergyBaseClass_Set(this,RegionMergingThresholdIn)

          IMPLICIT NONE

          CLASS(RCInternalEnergyBaseClass) :: this

          REAL(MK), INTENT(IN   ) :: RegionMergingThresholdIn
        END SUBROUTINE RCInternalEnergyBaseClass_Set
#endif

        FUNCTION DTYPE(RCInt_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel,e_merge)

          IMPLICIT NONE

          CLASS(RCInternalEnergyBaseClass)                     :: this

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER      :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER      :: imageIn
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),    POINTER       :: labelsIn
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),  POINTER       :: labelsIn
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

          DTYPE(RCInt_EvaluateEnergyDifference) =          &
          & this%DTYPE(EvaluateEnergyDifference)_(imageIn, &
          & labelsIn,coord,oldlabel,newlabel)

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
        &          imageIn,coord,aFromLabel,aToLabel,info)
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
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: imageIn
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: imageIn
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
          intensity=REAL(imageIn(coord(1),coord(2)),ppm_kind_double)
#elif __DIME == __3D
          intensity=REAL(imageIn(coord(1),coord(2),coord(3)),ppm_kind_double)
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

          LOGICAL :: sendrecv1,sendrecv2,lshrink

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

             nsize=this%size()
             ! Shrinkage ratio of 4 is a rule of thumb
             lshrink=nsize_.GT.4*nsize

             IF (lshrink) THEN
                CALL htable%destroy(info)
                or_fail("htable%destroy")

                CALL htable%create(nsize+1,info)
                or_fail("create htable")
             ENDIF

#ifdef __MPI
             CALL MPI_Wait(requestSums,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
#endif

             FORALL (i=0:nsize_,this%gCount(i).GT.halfd)
                this%gSums(i) =this%gSums(i) +this%lSumslSumsq(2*i)
                this%gSumsq(i)=this%gSumsq(i)+this%lSumslSumsq(2*i+1)
             END FORALL

             DEALLOCATE(this%lSumslSumsq,STAT=info)
             or_fail_dealloc("local array deallocation failed!")

             IF (lshrink) THEN
                CALL this%shrink(info)
                or_fail("this%shrink")

                DO i=0,nsize
                   CALL htable%insert(this%Rlabel(i),i,info)
                ENDDO
             ENDIF
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

          DO i=1,e_data%size()
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

          CALL htable%remove(label_region,info,.TRUE.)
          or_fail("hash remove")

          this%gCount(region_)=zerod
          this%gSums(region_) =zerod
          this%gSumsq(region_)=zerod
          this%Rlabel(region_)=-1

          end_subroutine()

        END SUBROUTINE RemoveFGRegion
#endif

        FUNCTION DTYPE(E_Gamma_EvaluateEnergyDifference)(this, &
        &        imageIn,labelsIn,coord,oldlabel,newlabel)

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          CLASS(E_Gamma)                                        :: this

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
          & labelsIn(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)
#elif __DIME == __3D
          & labelsIn(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)
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
          !!! Method is used to prepare energy functions.
          !!! It is called only once in the beginning of the filter.
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

          IF (ALLOCATED(this%NeighborsPoints)) THEN
             DEALLOCATE(this%NeighborsPoints,STAT=info)
             or_fail_dealloc("this%NeighborsPoints")
          ENDIF

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
        &        imageIn,labelsIn,coord,oldlabel,newlabel)

          IMPLICIT NONE

          CLASS(E_ContourLengthApprox)                          :: this

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
          tmplabels => labelsIn(coord(1)-e_lX:coord(1)+e_lX, &
          &                    coord(2)-e_lY:coord(2)+e_lY)
#elif __DIME == __3D
          tmplabels => labelsIn(coord(1)-e_lX:coord(1)+e_lX, &
          &                    coord(2)-e_lY:coord(2)+e_lY,  &
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

#if   __DIME == __3D
        ! Constructor
        SUBROUTINE ppm_rc_energy_parameter_redefine(info)
        !!! This subroutine will change all the energy variables
        !!! which are used in equilibrium phase

        IMPLICIT NONE
        INTEGER, INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_energy_parameter_redefine'
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        ! Now we set the running ghostsize to the correct ghost size
        ! after equilibrium
        ghostsize_run=ghostsize_equil

        DEALLOCATE(ghostsize_equil,STAT=info)
        or_fail_dealloc("ghostsize_equil")

        !-------------------------------------------------------------------------
        !  Initialize external energy related terms
        !-------------------------------------------------------------------------
        SELECT CASE (TRIM(energy_ext_name))
        CASE ("PC")
           CALL e_data%create(1,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")
        CASE ("PCGAUSSIAN")
           CALL e_data%create(2,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")
        CASE ("PCPOISSON")
           CALL e_data%create(3,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")
        CASE ("PS")
           CALL e_data%create(11,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           IF (energy_local_window_radius_equil.LT.bigs) THEN
              SELECT TYPE (e_data)
              TYPE IS (E_PS)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           ENDIF
        CASE ("PSGAUSSIAN")
           CALL e_data%create(12,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           IF (energy_local_window_radius_equil.LT.bigs) THEN
              SELECT TYPE (e_data)
              TYPE IS (E_PSGaussian)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           ENDIF
        CASE ("PSPOISSON")
           CALL e_data%create(13,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           IF (energy_local_window_radius_equil.LT.bigs) THEN
              SELECT TYPE (e_data)
              TYPE IS (E_PSPoisson)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           ENDIF
        END SELECT

        energy_coeff_balloon=energy_coeff_balloon_equil

        IF (energy_coeff_balloon.LT.bigs) THEN
           IF (energy_coeff_balloon.GT.smallest.OR.energy_coeff_balloon.LT.-smallest) THEN
              e_data%m_EnergyFunctional=e_data%m_EnergyFunctional+1000
           ENDIF
        ENDIF

        !-------------------------------------------------------------------------
        !  Initialize internal energy related terms
        !-------------------------------------------------------------------------
        SELECT CASE (TRIM(energy_int_name))
        CASE ("GAMMA")
           CALL e_length%create(31,energy_coeff_length,info)
           or_fail("e_length%create")
        CASE ("CURV")
           CALL e_length%create(32,energy_coeff_length,info)
           or_fail("e_length%create")

           IF (energy_curvature_mask_radius_equil.LT.bigs) THEN
              SELECT TYPE (e_length)
              TYPE IS (E_ContourLengthApprox)
                 e_length%m_Radius=energy_curvature_mask_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_length%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_length%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_length%PrepareEnergyCalculation")
              END SELECT
           ENDIF
        END SELECT

        energy_coeff_outward_flow=energy_coeff_outward_flow_equil

        IF (energy_coeff_outward_flow.LT.bigs) THEN
           IF (energy_coeff_outward_flow.GT.smallest.OR.energy_coeff_outward_flow.LT.-smallest) THEN
              e_length%m_EnergyFunctional=e_length%m_EnergyFunctional+1000
           ENDIF
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
        END SUBROUTINE ppm_rc_energy_parameter_redefine

        ! Constructor
        SUBROUTINE ppm_rc_energy_parameter_redefine_mcmc(info)
        !!! This subroutine will change all the energy variables
        !!! which are used in MCMC sampling

        IMPLICIT NONE
        INTEGER, INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_energy_parameter_redefine_mcmc'
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-------------------------------------------------------------------------
        ! Now we set the running ghostsize to the correct ghost size
        ! after equilibrium
        !-------------------------------------------------------------------------
        ghostsize_run=ghostsize_mcmc

        DEALLOCATE(ghostsize_mcmc,STAT=info)
        or_fail_dealloc("ghostsize_mcmc")

        !-------------------------------------------------------------------------
        !  If we are continuing after segmentation we migh need to update energy parameters
        !-------------------------------------------------------------------------
        IF (MCMCcontinue) THEN
           IF (energy_coeff_data_mcmc           .LT.bigs) energy_coeff_data           =energy_coeff_data_mcmc
           IF (energy_coeff_length_mcmc         .LT.bigs) energy_coeff_length         =energy_coeff_length_mcmc
           IF (energy_region_merge_ths_mcmc     .LT.bigs) energy_region_merge_ths     =energy_region_merge_ths_mcmc
           IF (energy_local_window_radius_mcmc  .LT.bigs) energy_local_window_radius  =energy_local_window_radius_mcmc
           IF (energy_curvature_mask_radius_mcmc.LT.bigs) energy_curvature_mask_radius=energy_curvature_mask_radius_mcmc
           IF (energy_coeff_balloon_mcmc        .LT.bigs) energy_coeff_balloon        =energy_coeff_balloon_mcmc
           IF (energy_coeff_outward_flow_mcmc   .LT.bigs) energy_coeff_outward_flow   =energy_coeff_outward_flow_mcmc

           !-------------------------------------------------------------------------
           !  Initialize external energy related terms
           !-------------------------------------------------------------------------
           SELECT CASE (TRIM(energy_ext_name))
           CASE ("PC")
              CALL e_data%create(1,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")
           CASE ("PCGAUSSIAN")
              CALL e_data%create(2,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")
           CASE ("PCPOISSON")
              CALL e_data%create(3,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")
           CASE ("PS")
              CALL e_data%create(11,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")

              SELECT TYPE (e_data)
              TYPE IS (E_PS)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           CASE ("PSGAUSSIAN")
              CALL e_data%create(12,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")

              SELECT TYPE (e_data)
              TYPE IS (E_PSGaussian)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           CASE ("PSPOISSON")
              CALL e_data%create(13,energy_coeff_data,info,energy_region_merge_ths)
              or_fail("e_data%create")

              SELECT TYPE (e_data)
              TYPE IS (E_PSPoisson)
                 e_data%m_Radius=energy_local_window_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_data%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_data%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_data%PrepareEnergyCalculation")
              END SELECT
           END SELECT

           IF (energy_coeff_balloon.LT.bigs) THEN
              IF (energy_coeff_balloon.GT.smallest.OR.energy_coeff_balloon.LT.-smallest) THEN
                 e_data%m_EnergyFunctional=e_data%m_EnergyFunctional+1000
              ENDIF
           ENDIF

           !-------------------------------------------------------------------------
           !  Initialize internal energy related terms
           !-------------------------------------------------------------------------
           SELECT CASE (TRIM(energy_int_name))
           CASE ("GAMMA")
              CALL e_length%create(31,energy_coeff_length,info)
              or_fail("e_length%create")
           CASE ("CURV")
              CALL e_length%create(32,energy_coeff_length,info)
              or_fail("e_length%create")

              SELECT TYPE (e_length)
              TYPE IS (E_ContourLengthApprox)
                 e_length%m_Radius=energy_curvature_mask_radius

                 SELECT CASE (ppm_rc_dim)
                 CASE (2)
                    CALL e_length%PrepareEnergyCalculation_2d(info)
                 CASE (3)
                    CALL e_length%PrepareEnergyCalculation_3d(info)
                 END SELECT
                 or_fail("e_length%PrepareEnergyCalculation")
              END SELECT
           END SELECT

           IF (energy_coeff_outward_flow.LT.bigs) THEN
              IF (energy_coeff_outward_flow.GT.smallest.OR.energy_coeff_outward_flow.LT.-smallest) THEN
                 e_length%m_EnergyFunctional=e_length%m_EnergyFunctional+1000
              ENDIF
           ENDIF
        ENDIF !(MCMCcontinue)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
        END SUBROUTINE ppm_rc_energy_parameter_redefine_mcmc
#endif