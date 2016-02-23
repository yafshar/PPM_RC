        ! Set image values to a user-specified value if they are below,
        ! above, or between simple threshold values.
        !
        ! ThresholdImageFilter sets image values to a user-specified "outside"
        ! value (by default, "black") if the image values are below, above, or
        ! between simple threshold values.
        !
        ! You can use the MaskValue to convert the threshold image into a Mask.
        TYPE :: ThresholdImageFilter
          REAL(ppm_kind_double) :: m_OutsideValue
          REAL(ppm_kind_double) :: m_Lower
          REAL(ppm_kind_double) :: m_Upper
        CONTAINS
          PROCEDURE :: ThresholdAbove   => ThresholdImageFilter_ThresholdAbove
          !!!Set the threshold parameters so it will
          !!!threshold above thresh value
          PROCEDURE :: ThresholdBelow   => ThresholdImageFilter_ThresholdBelow
          !!!Set the threshold parameters so it will
          !!!threshold below thresh value
          PROCEDURE :: ThresholdOutside => ThresholdImageFilter_ThresholdOutside
          !!!Set the threshold parameters so it will
          !!!threshold below and above thresh value
          PROCEDURE :: GenerateData_2d  => ThresholdImageFilter_GenerateData_2d
          PROCEDURE :: GenerateData_3d  => ThresholdImageFilter_GenerateData_3d
          !!!generate the threshold image, using the MaskValue you can
          !!!produce a Mask from threshold image
        END TYPE ThresholdImageFilter

        !Histogram
        !This class stores measurement vectors in the context of histogram.
        !
        !Histogram bins can be regularly spaced. The storage for the histogram is
        !managed via the FrequencyContainer & OffsetTable.
        !
        !Frequencies of a bin index.
        !
        !Measurements can be queried by bin index.
        !In this case, the measurement returned in the centroid
        !of the histogram bin.
        !
        !The Initialize() method is used to specified the number of bins for
        !the histogram for regularly spaced bins to defined.
        !
        TYPE,ABSTRACT :: Histogram
          INTEGER                                          :: histSize
          !!!The number of bins

          REAL(ppm_kind_double)                            :: histOffset
          !!!The interval of each bin
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: FrequencyContainer
          !!!Frequencies of a bin index.
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: OffsetTable

        CONTAINS
          PROCEDURE :: Initialize => Histogram_Initialize
          PROCEDURE :: create_2d  => Histogram_create_2d
          PROCEDURE :: create_3d  => Histogram_create_3d
          PROCEDURE :: destroy    => Histogram_destroy
          PROCEDURE :: Mean       => Histogram_Mean
          !!!Get the mean value
          PROCEDURE :: Quantile   => Histogram_Quantile
          !!!Get the pth percentile value.
          !!!
          !!!Let assume n = the index of the bin where the p-th percentile value is,
          !!!interval = histOffset ,
          !!!pp = cumlated proportion until n-1 bin;
          !!!and pb = frequency of the bin / total frequency of the dimension.
          !!!
          !!!If p is less than 0.5,
          !!!the percentile value =
          !!!min + ((p - pp ) / pb) * interval
          !!!If p is greater than or equal to 0.5
          !!!the percentile value =
          !!!max - ((pp - p) / pb) * interval  */

        END TYPE Histogram


        !HistogramThresholdCalculator
        !Base class to compute a threshold value based on the histogram of an image
        TYPE,EXTENDS(Histogram),ABSTRACT :: HistogramThresholdCalculator
        END TYPE HistogramThresholdCalculator

        !HistogramThresholdImageFilter
        !Threshold an image using a HistogramThresholdCalculator
        !
        !This filter creates a binary thresholded image that separates an
        !image into foreground and background components. The filter
        !computes the threshold using a user provided HistogramThresholdCalculator and
        !applies that theshold to the input image using the
        !BinaryThresholdImageFilter.
        !
        TYPE,ABSTRACT :: HistogramThresholdImageFilter
          TYPE(ThresholdImageFilter) :: Calculator
        END TYPE HistogramThresholdImageFilter

        !OtsuThresholdCalculator
        !brief Computes the Otsu's threshold for an image.
        !
        !This calculator computes the Otsu's threshold which separates an image
        !into foreground and background components. The method relies on a
        !histogram of image intensities. The basic idea is to maximize the
        !between-class variance.
        TYPE,EXTENDS(HistogramThresholdCalculator) :: OtsuThresholdCalculator
        CONTAINS
          PROCEDURE :: SetInput_2d => OtsuThresholdCalculator_SetInput_2d
          PROCEDURE :: SetInput_3d => OtsuThresholdCalculator_SetInput_3d
        END TYPE OtsuThresholdCalculator


        !OtsuThresholdImageFilter
        !Threshold an image using the Otsu Threshold
        !
        !This filter creates a binary thresholded image that separates an
        !image into foreground and background components. The filter
        !computes the threshold using the OtsuThresholdCalculator and
        !applies that theshold to the input image using the
        !BinaryThresholdImageFilter.
        TYPE,EXTENDS(HistogramThresholdImageFilter) :: OtsuThresholdImageFilter
        CONTAINS
          PROCEDURE :: GenerateData_2d => OtsuThresholdImageFilter_GenerateData_2d
          PROCEDURE :: GenerateData_3d => OtsuThresholdImageFilter_GenerateData_3d
        END TYPE OtsuThresholdImageFilter

        !OtsuMultipleThresholdsCalculator
        !Computes Otsu's thresholds for a histogram.
        !
        !You plug in the target histogram using SetInputHistogram method and
        !specify the number of thresholds you want to be computed. Then call
        !the GenerateData method to run the alogithm.
        !
        !The thresholds are computed so that the between-class variance is
        !maximized.

        !OtsuMultipleThresholdsImageFilter
        !Threshold an image using multiple Otsu Thresholds.
        !
        !
        !This filter creates a labeled image that separates the input
        !image into various classes. The filter computes the
        !thresholds using the OtsuMultipleThresholdsCalculator and
        !applies those thesholds to the input image.
        !
        TYPE,EXTENDS(HistogramThresholdImageFilter) :: OtsuMultipleThresholdsImageFilter
        CONTAINS
          PROCEDURE :: GenerateData_2d => OtsuMultipleThresholdsImageFilter_GenerateData_2d
          PROCEDURE :: GenerateData_3d => OtsuMultipleThresholdsImageFilter_GenerateData_3d
        END TYPE OtsuMultipleThresholdsImageFilter


        ! RecursiveSeparableImageFilter is the base class for recursive
        ! filters that are applied in each dimension separately.
        TYPE,ABSTRACT :: RecursiveSeparableImageFilter
          ! Causal coefficients that multiply the input data.
          REAL(MK) :: N0
          REAL(MK) :: N1
          REAL(MK) :: N2
          REAL(MK) :: N3
          ! Recursive coefficients that multiply previously computed values
          ! at the output. These are the same for the causal and
          ! anti-causal parts of the filter.
          REAL(MK) :: D1
          REAL(MK) :: D2
          REAL(MK) :: D3
          REAL(MK) :: D4
          ! Anti-causal coefficients that multiply the input data. */
          REAL(MK) :: M1
          REAL(MK) :: M2
          REAL(MK) :: M3
          REAL(MK) :: M4
          ! Recursive coefficients to be used at the boundaries to simulate
          ! edge extension boundary conditions. */
          REAL(MK) :: BN1
          REAL(MK) :: BN2
          REAL(MK) :: BN3
          REAL(MK) :: BN4

          REAL(MK) :: BM1
          REAL(MK) :: BM2
          REAL(MK) :: BM3
          REAL(MK) :: BM4

          LOGICAL  :: NormalizeAcrossScale=.FALSE.
        CONTAINS
          PROCEDURE :: FilterDataArray_2d
          PROCEDURE :: FilterDataArray_3d
          GENERIC   :: FilterDataArray =>  &
          &            FilterDataArray_2d, &
          &            FilterDataArray_3d
        END TYPE RecursiveSeparableImageFilter

        ! RecursiveGaussianImageFilter is the base class for recursive filters that
        ! approximate convolution with the Gaussian kernel.
        !
        ! This class implements the recursive filtering
        ! method proposed by R.Deriche in IEEE-PAMI
        ! Vol.12, No.1, January 1990, pp 78-87,
        ! "Fast Algorithms for Low-Level Vision"
        !
        ! Details of the implementation are described in the technical report:
        ! R. Deriche, "Recursively Implementing The Gaussian and Its Derivatives",
        ! INRIA, 1993, ftp://ftp.inria.fr/INRIA/tech-reports/RR/RR-1893.ps.gz
        !
        TYPE,EXTENDS(RecursiveSeparableImageFilter) :: RecursiveGaussianImageFilter
          REAL(MK) :: Sigma=-one
          !!!Sigma of the gaussian kernel
          INTEGER  :: Order=-1
          !!!zeroth, first or second order derivative.
          !!!ZeroOrder   --- 0
          !!!FirstOrder  --- 1
          !!!SecondOrder --- 2
        CONTAINS
          PROCEDURE :: ComputeNCoefficients
          PROCEDURE :: ComputeDCoefficients
          PROCEDURE :: ComputeRemainingCoefficients
          PROCEDURE :: SetUp => RecursiveGaussianImageFilter_SetUp
        END TYPE RecursiveGaussianImageFilter




























