
















BlobDetectionImageFilter

LaplacianRecursiveGaussianImageFilter :: LaplacianFilterType
ScalarImageToHistogramGenerator :: HistogramGeneratorType
HistogramGeneratorType::HistogramType HistogramType
OtsuMultipleThresholdsCalculator :: OtsuThsCalcType
MinimumProjectionImageFilter :: MinProjIFType
ThresholdImageFilter :: ThresholdImageFilterType
FlatStructuringElement :: SEType
RegionalMinimaImageFilter :: RegMinImageFilterType
BinaryDilateImageFilter :: DilateImageFilterType   m_DilateFilter
AndImageFilter :: AndFilterType
AddImageFilter :: AddFilterType
BinaryThresholdImageFilter :: BinaryThsFilterType





  unsigned int m_NumberOfScales;
  float m_MaxDiameter;
  float m_MinDiameter;



  BlobDetectionImageFilter<TInputImage, TOutputImage>::BlobDetectionImageFilter()
{
    m_MaxDiameter = 20;
    m_MinDiameter = 5;
    m_NumberOfScales = 5;

    m_DilateFilter = DilateImageFilterType::New();
}

BlobDetectionImageFilter :: GenerateData
{
  ! Allocate the scale space image (DIM+1)
  RealCoDimImageType::Pointer vScaleSpace = RealCoDimImageType::New();
  RealCoDimImageType::SizeType vSize;
  RealCoDimImageType::RegionType vRegion;
  RealImageType::RegionType vSliceRegion = this->GetInput()->GetRequestedRegion();
  for(unsigned int vD =0; vD < TInputImage::ImageDimension; vD++) {
    vSize[vD] = vSliceRegion.GetSize()[vD];
  }
  vSize[TInputImage::ImageDimension] = m_NumberOfScales; ! set the last
  vRegion.SetSize(vSize);
  vScaleSpace->SetRegions(vRegion);
  vScaleSpace->Allocate();

  LaplacianFilterType:: vLaplacianFilter = LaplacianFilterType::New();
  vLaplacianFilter->SetInput(this->GetInput());
  vLaplacianFilter->SetNormalizeAcrossScale(true);

  itk::ImageRegionIterator<RealCoDimImageType> vCopyIt(vScaleSpace, vScaleSpace->GetLargestPossibleRegion());


    ! Generate the scale space image
  float vSigmaMax = 0.5f*m_MaxDiameter/1.414f;
  float vSigmaMin = 0.5f*m_MinDiameter/1.414f;
  float vSigmaInc = (vSigmaMax-vSigmaMin)/(m_NumberOfScales-1);
  float vSigma = vSigmaMin;
  for(unsigned int vI = 0; vI < m_NumberOfScales; vI++) {
    std::cout << "Sigma:\t" <<  vSigma << std::endl;

    vLaplacianFilter->SetSigma(vSigma);
    vLaplacianFilter->Update();

    ! copy the slice to the scale space image
    itk::ImageRegionConstIterator<RealImageType> vOrigIt(vLaplacianFilter->GetOutput(), vLaplacianFilter->GetOutput()->GetBufferedRegion());

    for(vOrigIt.GoToBegin(); !vOrigIt.IsAtEnd(); ++vOrigIt) {
      vCopyIt.Set(vOrigIt.Get());
      ++vCopyIt;
    }
    vSigma += vSigmaInc;

  }

  ! First do a projection, its faster and the thresholding does not
  ! depend on the number of scale stages.
  MinProjIFType:: vMinProjFilter = MinProjIFType::New();
  vMinProjFilter->SetInput(vScaleSpace);
  vMinProjFilter->Update();

  !Find the threshold to select only meaningful/strong local minima
  HistogramGeneratorType:: vHistogramGenerator = HistogramGeneratorType::New();
  vHistogramGenerator->SetInput(vMinProjFilter->GetOutput());
  unsigned int vNumberOfBins = 256;
  vHistogramGenerator->SetNumberOfBins(vNumberOfBins);
  vHistogramGenerator->Compute();

  OtsuThsCalcType:: vOtsuThsCalculator = OtsuThsCalcType::New();
  vOtsuThsCalculator->SetInputHistogram(vHistogramGenerator->GetOutput());
  unsigned int vNumberOfThresholds = 2;
  vOtsuThsCalculator->SetNumberOfThresholds(vNumberOfThresholds);
  vOtsuThsCalculator->Update();
  OtsuThsCalcType::OutputType vThs = vOtsuThsCalculator->GetOutput();

  ThresholdImageFilterType:: vThsIF = ThresholdImageFilterType::New();
  vThsIF->SetInput(vMinProjFilter->GetOutput());
  vThsIF->SetUpper(vThs[vNumberOfThresholds-1]);

  ! Within the thresholded image, find local minima
  RegMinImageFilterType:: vRegMinFilter = RegMinImageFilterType::New();
  vRegMinFilter->SetInput(vThsIF->GetOutput());
  vRegMinFilter->SetFullyConnected(true);
  vRegMinFilter->SetForegroundValue(1);

  ! Create a small blob for each local minimum.
  m_DilateFilter->SetDilateValue(1);

  SEType::RadiusType vSERadius;
  ! TODO: the radius of each blob can be found in the scale space image.
  vSERadius.Fill(m_MinDiameter/3.0f); ! such that discs of different blobs do not touch.
  m_DilateFilter->SetKernel(SEType::Ball(vSERadius));
  m_DilateFilter->SetInput(vRegMinFilter->GetOutput());

  ! graft the output.
  m_DilateFilter->GraftOutput( this->GetOutput() );
  m_DilateFilter->Update();
  this->GraftOutput(m_DilateFilter->GetOutput() );
}












