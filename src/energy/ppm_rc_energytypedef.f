        ! Base class for energy functions for
        ! the Region Competition optimizer and sampler

        TYPE, ABSTRACT :: RCEnergyBaseClass
          REAL(ppm_kind_double) :: m_Coefficient
          !!!Energy term coefficient

          ! 1 ---> E_PC
          ! 2 ---> E_PCGaussian
          ! 3 ---> E_PCPoisson
          !
          !11 ---> E_PS
          !12 ---> E_PSGaussian
          !13 ---> E_PSPoisson
          !
          !31 ---> E_Gamma
          !32 ---> E_ContourLengthApprox
          !
          INTEGER  :: m_EnergyFunctional
        CONTAINS
          PROCEDURE :: create => RCEnergyBaseClass_create

          PROCEDURE(RCEnergyBaseClass_Set), DEFERRED :: Set

          PROCEDURE(RC_EvaluateEnergyDifference_2d), &
          &          DEFERRED :: EvaluateEnergyDifference_2d
          PROCEDURE(RC_EvaluateEnergyDifference_3d), &
          &          DEFERRED :: EvaluateEnergyDifference_3d
          GENERIC :: EvaluateEnergyDifference =>  &
          &          EvaluateEnergyDifference_2d, &
          &          EvaluateEnergyDifference_3d

          PROCEDURE :: CalculateVariance
          PROCEDURE :: CalculateScaledSphereVolume
        END TYPE RCEnergyBaseClass

        TYPE, EXTENDS(RCEnergyBaseClass), ABSTRACT :: RCExternalEnergyBaseClass
          REAL(MK)                                         :: RegionMergingThreshold

          !----------------------------------------------------------------------
          !  Region statistics.
          !  This is global knowledge, replicated on all processors.
          !----------------------------------------------------------------------
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: gCount
          !!! The number of pixels within a region
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: gSums
          !!! Sum of intensities within a region - used to compute means
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: gSumsq
          !!! Sum of squares of intensities within regions - used to compute variances
          INTEGER,               DIMENSION(:), ALLOCATABLE :: Rlabel

          !----------------------------------------------------------------------
          !  Region statistics.
          !  This is local knowledge. Only on local processor.
          !----------------------------------------------------------------------
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: lCount
          !!! The number of pixels within a local regin
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: lSumslSumsq
          !!! Sum of intensities within a local region - used to compute means
          !!! Sum of squares of intensities within local regions - used to compute variances
        CONTAINS
          PROCEDURE :: Set => RCExternalEnergyBaseClass_Set

          PROCEDURE :: lalloc  => lalloc_RCExternalEnergyBaseClass
          PROCEDURE :: galloc  => galloc_RCExternalEnergyBaseClass
          PROCEDURE :: grow    => grow_RCExternalEnergyBaseClass
          PROCEDURE :: destroy => destroy_RCExternalEnergyBaseClass

          PROCEDURE(destroy_RCExternalEnergyBaseClass_), DEFERRED :: destroy_

          PROCEDURE :: EvaluateEnergyDifference_2d => &
          &            RCExt_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d => &
          &            RCExt_EvaluateEnergyDifference_3d

          PROCEDURE(RCExt_EvaluateEnergyDifference_2d_), &
          &            DEFERRED :: EvaluateEnergyDifference_2d_
          PROCEDURE(RCExt_EvaluateEnergyDifference_3d_), &
          &            DEFERRED :: EvaluateEnergyDifference_3d_

          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            RCExt_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            RCExt_EvaluateEnergyDifference_E_Merge_3d

          PROCEDURE :: CalculateTotalEnergy_2d
          PROCEDURE :: CalculateTotalEnergy_3d
          GENERIC   :: CalculateTotalEnergy => &
          &            CalculateTotalEnergy_2d, &
          &            CalculateTotalEnergy_3d

          !----------------------------------------------------------------------
          !  Kullback–Leibler divergence is a non-symmetric measure of the
          !  difference between two probability distributions.
          !----------------------------------------------------------------------
          PROCEDURE :: CalculateKullbackLeiblerDistance

          PROCEDURE :: AddPoint_2d
          PROCEDURE :: AddPoint_3d
          PROCEDURE :: AddPoint_i
          PROCEDURE :: AddPoint_rs
          PROCEDURE :: AddPoint_r
          GENERIC   :: AddPoint =>  &
          &            AddPoint_2d, &
          &            AddPoint_3d, &
          &            AddPoint_i,  &
          &            AddPoint_rs, &
          &            AddPoint_r

          PROCEDURE :: RemovePoint_2d
          PROCEDURE :: RemovePoint_3d
          PROCEDURE :: RemovePoint_i
          PROCEDURE :: RemovePoint_rs
          PROCEDURE :: RemovePoint_r
          GENERIC   :: RemovePoint =>  &
          &            RemovePoint_2d, &
          &            RemovePoint_3d, &
          &            RemovePoint_i,  &
          &            RemovePoint_rs, &
          &            RemovePoint_r

          PROCEDURE :: KillRegion

          PROCEDURE :: UpdateStatisticsWhenJump_2d
          PROCEDURE :: UpdateStatisticsWhenJump_3d
          GENERIC   :: UpdateStatisticsWhenJump =>  &
          &            UpdateStatisticsWhenJump_2d, &
          &            UpdateStatisticsWhenJump_3d

          PROCEDURE :: UpdateStatistics

          PROCEDURE :: RemoveNotSignificantRegions
          procedure :: RemoveFGRegion
        END TYPE RCExternalEnergyBaseClass

        TYPE, EXTENDS(RCEnergyBaseClass), ABSTRACT :: RCInternalEnergyBaseClass
        CONTAINS
          PROCEDURE :: Set => RCInternalEnergyBaseClass_Set

          PROCEDURE :: EvaluateEnergyDifference_2d => &
          &            RCInt_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d => &
          &            RCInt_EvaluateEnergyDifference_3d

          PROCEDURE(RCInt_EvaluateEnergyDifference_2d_), &
          &            DEFERRED :: EvaluateEnergyDifference_2d_
          PROCEDURE(RCInt_EvaluateEnergyDifference_3d_), &
          &            DEFERRED :: EvaluateEnergyDifference_3d_
        END TYPE RCInternalEnergyBaseClass

        !----------------------------------------------------------------------
        !  Energy terms which will be used in the segmentation
        !----------------------------------------------------------------------

        TYPE, EXTENDS(RCInternalEnergyBaseClass) :: E_Gamma
        CONTAINS
          PROCEDURE :: EvaluateEnergyDifference_2d_ => &
          &            E_Gamma_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ => &
          &            E_Gamma_EvaluateEnergyDifference_3d
        END TYPE E_Gamma

        TYPE, EXTENDS(RCInternalEnergyBaseClass) :: E_ContourLengthApprox
        !!! Returns A regularizing energy difference based
        !!! on local curvature computation.

        !!! For details, see "J. Kybic and J. Kratky, "Discrete curvature calculation
        !!! for fast level set segmentation," in Proc. 16th IEEE Int. Conf. Image
        !!! Process., Nov. 2009, pp. 3017-3020
          REAL(MK)                             :: m_Radius
          REAL(MK)                             :: m_Volume
          REAL(MK)                             :: m_Prefactor

          INTEGER                              :: NumberOfNeighbors
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: NeighborsPoints
        CONTAINS
          PROCEDURE :: destroy => E_ContourLengthApprox_destroy

          PROCEDURE :: PrepareEnergyCalculation_2d =>  &
          &            E_ContourLengthApprox_PrepareEnergyCalculation_2d
          PROCEDURE :: PrepareEnergyCalculation_3d =>  &
          &            E_ContourLengthApprox_PrepareEnergyCalculation_3d

          PROCEDURE :: EvaluateEnergyDifference_2d_ => &
          &            E_ContourLengthApprox_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ => &
          &            E_ContourLengthApprox_EvaluateEnergyDifference_3d
        END TYPE E_ContourLengthApprox

        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PC
        !!! Computes energy differences for a pixel change in for a piece-wise
        !!! constant image model.
        !!!
        !!! The energy to minimize is \f$E = \sum_i^M(\mu_i - I(x))^2\f$ with M being the
        !!! number of regions, I the image and \f$\mu_i\f$ the mean of region i.
        !!!
        !!! In case of 2 regions this image model is also called Chan-Vese model.
        CONTAINS
          PROCEDURE :: destroy_ => E_PC_destroy

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PC_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PC_EvaluateEnergyDifference_3d

          PROCEDURE :: CalculateTotalEnergy_2d => &
          &            E_PC_CalculateTotalEnergy_2d
          PROCEDURE :: CalculateTotalEnergy_3d => &
          &            E_PC_CalculateTotalEnergy_3d

        END TYPE E_PC

        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PCGaussian
        !!! Computes energy differences for a pixel change in for a piece-wise
        !!! constant image model and with i.i.d. Gaussian noise.
        !!!
        !!! The energy to minimize is \f$E = \sum_i^M(\mu_i - I(x))^2\f$ with M being the
        !!! number of regions, I the image and \f$\mu_i\f$ the mean of region i.
        !!!
        !!! In case of 2 regions this image model is also called Chan-Vese model.
        CONTAINS
          PROCEDURE :: destroy_ => E_PCGaussian_destroy

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PCGaussian_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PCGaussian_EvaluateEnergyDifference_3d
          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            E_PCGaussian_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            E_PCGaussian_EvaluateEnergyDifference_E_Merge_3d

        END TYPE E_PCGaussian

        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PCPoisson
        !!! Computes energy differences for a pixel change in for a piece-wise
        !!! constant image model and with i.i.d. Poisson noise.
        !!!
        !!! The energy to minimize is \f$E = \sum_i^M(\mu_i - I(x))^2\f$ with M being the
        !!! number of regions, I the image and \f$\mu_i\f$ the mean of region i.
        !!!
        !!! In case of 2 regions this image model is also called Chan-Vese model.
        CONTAINS
          PROCEDURE :: destroy_ => E_PCPoisson_destroy

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PCPoisson_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PCPoisson_EvaluateEnergyDifference_3d
          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            E_PCPoisson_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            E_PCPoisson_EvaluateEnergyDifference_E_Merge_3d

        END TYPE E_PCPoisson


        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PS
        !!! Computes energy differences for a pixel change in
        !!! for a piece-wise smooth image model.
        !!!
        !!! The energy to minimize is ...
          REAL(MK)                             :: m_Radius

          INTEGER                              :: NumberOfNeighbors
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: NeighborsPoints

        CONTAINS
          PROCEDURE :: destroy_ => E_PS_destroy

          PROCEDURE :: PrepareEnergyCalculation_2d => &
          &            E_PS_PrepareEnergyCalculation_2d
          PROCEDURE :: PrepareEnergyCalculation_3d => &
          &            E_PS_PrepareEnergyCalculation_3d

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PS_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PS_EvaluateEnergyDifference_3d

          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            E_PS_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            E_PS_EvaluateEnergyDifference_E_Merge_3d

          PROCEDURE :: CalculateTotalEnergy_2d => &
          &            E_PS_CalculateTotalEnergy_2d
          PROCEDURE :: CalculateTotalEnergy_3d => &
          &            E_PS_CalculateTotalEnergy_3d

        END TYPE E_PS

        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PSGaussian
        !!! Computes energy differences for a pixel change in for a piece-wise
        !!! smooth image model and with i.i.d. Gaussian noise.
        !!!
        !!! The energy to minimize is ...
          REAL(MK)                             :: m_Radius

          INTEGER                              :: NumberOfNeighbors
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: NeighborsPoints

        CONTAINS
          PROCEDURE :: destroy_ => E_PSGaussian_destroy

          PROCEDURE :: PrepareEnergyCalculation_2d => &
          &            E_PSGaussian_PrepareEnergyCalculation_2d
          PROCEDURE :: PrepareEnergyCalculation_3d => &
          &            E_PSGaussian_PrepareEnergyCalculation_3d

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PSGaussian_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PSGaussian_EvaluateEnergyDifference_3d

          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            E_PSGaussian_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            E_PSGaussian_EvaluateEnergyDifference_E_Merge_3d
        END TYPE E_PSGaussian

        TYPE, EXTENDS(RCExternalEnergyBaseClass) :: E_PSPoisson
        !!! Computes energy differences for a pixel change in for a piece-wise
        !!! smooth image model and with i.i.d. Poisson noise.
        !!!
        !!! The energy to minimize is ...
          REAL(MK)                             :: m_Radius

          INTEGER                              :: NumberOfNeighbors
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: NeighborsPoints

        CONTAINS
          PROCEDURE :: destroy_ => E_PSPoisson_destroy

          PROCEDURE :: PrepareEnergyCalculation_2d => &
          &            E_PSPoisson_PrepareEnergyCalculation_2d
          PROCEDURE :: PrepareEnergyCalculation_3d => &
          &            E_PSPoisson_PrepareEnergyCalculation_3d

          PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
          &            E_PSPoisson_EvaluateEnergyDifference_2d
          PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
          &            E_PSPoisson_EvaluateEnergyDifference_3d

          PROCEDURE :: EvaluateEnergyDifference_2d__ => &
          &            E_PSPoisson_EvaluateEnergyDifference_E_Merge_2d
          PROCEDURE :: EvaluateEnergyDifference_3d__ => &
          &            E_PSPoisson_EvaluateEnergyDifference_E_Merge_3d
        END TYPE E_PSPoisson


    !    TYPE, EXTENDS(RCExternalEnergyBaseClass), ABSTRACT :: PCDeconvolutionEnergyBaseClass
    !    !!! Computes energy differences for a pixel change in for a piece-wise
    !    !!! constant deconvolution image model.

    !    InternalImagePointerType m_DeconvolutionModelImage
    !    InternalImagePointerType m_PSF
    !    DataImageConstPointerType m_PSFInput
    !    LOGICAL :: m_DoNormalizePSF = .FALSE.
    !    LOGICAL :: m_SwitchPointMode = .FALSE.
    !    LOGICAL :: m_OptimizationMode = .FALSE.


    !    CONTAINS
    !      PROCEDURE :: destroy_ => E_PC_destroy

    !      PROCEDURE :: EvaluateEnergyDifference_2d_ =>  &
    !      &            E_PC_EvaluateEnergyDifference_2d
    !      PROCEDURE :: EvaluateEnergyDifference_3d_ =>  &
    !      &            E_PC_EvaluateEnergyDifference_3d
    !    END TYPE PCDeconvolutionEnergyBaseClass



        !----------------------------------------------------------------------
        !  INTERFACES
        !----------------------------------------------------------------------
        INTERFACE
          ! Constructor
          SUBROUTINE RCEnergyBaseClass_Set(this,RegionMergingThreshold_)
            IMPORT                   :: RCEnergyBaseClass
            IMPORT                   :: MK
            IMPLICIT NONE
            CLASS(RCEnergyBaseClass) :: this
            REAL(MK),  INTENT(IN   ) :: RegionMergingThreshold_
          END SUBROUTINE RCEnergyBaseClass_Set

          FUNCTION RC_EvaluateEnergyDifference_2d(this, &
          &        image_,labels_,coord,oldlabel,newlabel,e_merge)
            IMPORT                                              :: RCEnergyBaseClass
            IMPORT                                              :: MK
            IMPLICIT NONE
            CLASS(RCEnergyBaseClass)                            :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
            INTEGER,                              INTENT(IN   ) :: oldlabel
            INTEGER,                              INTENT(IN   ) :: newlabel
            LOGICAL,                              INTENT(INOUT) :: e_merge
            REAL(MK)                                            :: RC_EvaluateEnergyDifference_2d
          END FUNCTION RC_EvaluateEnergyDifference_2d

          FUNCTION RC_EvaluateEnergyDifference_3d(this, &
          &        image_,labels_,coord,oldlabel,newlabel,e_merge)
            IMPORT                                                :: RCEnergyBaseClass
            IMPORT                                                :: MK
            IMPLICIT NONE
            CLASS(RCEnergyBaseClass)                              :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
            INTEGER,                                INTENT(IN   ) :: oldlabel
            INTEGER,                                INTENT(IN   ) :: newlabel
            LOGICAL,                                INTENT(INOUT) :: e_merge
            REAL(MK)                               :: RC_EvaluateEnergyDifference_3d
          END FUNCTION RC_EvaluateEnergyDifference_3d

          FUNCTION RCExt_EvaluateEnergyDifference_2d_(  &
          &        this,image_,labels_,coord,oldlabel,  &
          &        newlabel,oldlabel_region,newlabel_region)
            IMPORT                                              :: RCExternalEnergyBaseClass
            IMPORT                                              :: MK
            IMPLICIT NONE
            CLASS(RCExternalEnergyBaseClass)                    :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
            INTEGER,                              INTENT(IN   ) :: oldlabel
            INTEGER,                              INTENT(IN   ) :: newlabel
            INTEGER,                              INTENT(IN   ) :: oldlabel_region
            INTEGER,                              INTENT(IN   ) :: newlabel_region
            REAL(MK)                                            :: RCExt_EvaluateEnergyDifference_2d_
          END FUNCTION RCExt_EvaluateEnergyDifference_2d_

          FUNCTION RCExt_EvaluateEnergyDifference_3d_( &
          &        this,image_,labels_,coord,oldlabel, &
          &        newlabel,oldlabel_region,newlabel_region)
            IMPORT                                                :: RCExternalEnergyBaseClass
            IMPORT                                                :: MK
            IMPLICIT NONE
            CLASS(RCExternalEnergyBaseClass)                      :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
            INTEGER,                                INTENT(IN   ) :: oldlabel
            INTEGER,                                INTENT(IN   ) :: newlabel
            INTEGER,                                INTENT(IN   ) :: oldlabel_region
            INTEGER,                                INTENT(IN   ) :: newlabel_region
            REAL(MK)                                              :: RCExt_EvaluateEnergyDifference_3d_
          END FUNCTION RCExt_EvaluateEnergyDifference_3d_

          SUBROUTINE destroy_RCExternalEnergyBaseClass_(this,info)
            IMPORT                           :: RCExternalEnergyBaseClass
            IMPLICIT NONE
            CLASS(RCExternalEnergyBaseClass) :: this
            INTEGER,           INTENT(  OUT) :: info
          END SUBROUTINE destroy_RCExternalEnergyBaseClass_

          FUNCTION RCInt_EvaluateEnergyDifference_2d_(this, &
          &        image_,labels_,coord,oldlabel,newlabel)
            IMPORT                                              :: RCInternalEnergyBaseClass
            IMPORT                                              :: MK
            IMPLICIT NONE
            CLASS(RCInternalEnergyBaseClass)                    :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
            INTEGER,                              INTENT(IN   ) :: oldlabel
            INTEGER,                              INTENT(IN   ) :: newlabel
            REAL(MK)                               :: RCInt_EvaluateEnergyDifference_2d_
          END FUNCTION RCInt_EvaluateEnergyDifference_2d_

          FUNCTION RCInt_EvaluateEnergyDifference_3d_(this, &
          &        image_,labels_,coord,oldlabel,newlabel)
            IMPORT                                                :: RCInternalEnergyBaseClass
            IMPORT                                                :: MK
            IMPLICIT NONE
            CLASS(RCInternalEnergyBaseClass)                      :: this
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
            INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
            INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
            INTEGER,                                INTENT(IN   ) :: oldlabel
            INTEGER,                                INTENT(IN   ) :: newlabel
            REAL(MK)                                              :: RCInt_EvaluateEnergyDifference_3d_
          END FUNCTION RCInt_EvaluateEnergyDifference_3d_
        END INTERFACE
