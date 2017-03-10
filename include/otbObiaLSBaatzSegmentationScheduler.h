#ifndef otbObiaLSBaatzSegmentationScheduler_h
#define otbObiaLSBaatzSegmentationScheduler_h
#include "otbObiaBaatzSegmentationFilter.h"
#include "otbObiaLSImageToGraphScheduler.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaGenericRegionMergingFilter.h"

namespace otb
{
namespace obia
{

template< class TInputImage >
  class LSBaatzSegmentationScheduler : public LSImageToGraphScheduler< TInputImage,
                        typename ImageToBaatzGraphFilter<TInputImage>::OutputGraphType >
{
public:

  /** Standard class alias */
  using Self 		 = LSBaatzSegmentationScheduler;
  using SuperClass   = itk::Object;
  using Pointer      = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >; 

  /** Some convenient alias */
  using InputImageType                   = TInputImage;
  using InputInternalPixelType           = typename InputImageType::InternalPixelType;
  using InputImageReaderType             = ImageFileReader<InputImageType>;
  using MultiChannelExtractROIFilterType = MultiChannelExtractROI<InputInternalPixelType,InputInternalPixelType>;
  using ImageToBaatzGraphFilterType 	   = ImageToBaatzGraphFilter<InputImageType>;
  using GraphType                        = typename ImageToBaatzGraphFilterType::OutputGraphType;
  using NodeType                         = typename GraphType::NodeType;
  using GraphPointerType                 = typename GraphType::Pointer;
  using GraphOperationsType              = GraphOperations<GraphType>;
  using BaatzSegmentationFilterType	     = GenericRegionMergingFilter<GraphType, GraphType,
		  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  BaatzMergingCost<float, GraphType> ,
																	  BaatzHeuristic<GraphType>,
																	  BaatzUpdateAttribute<GraphType> >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LSBaatzSegmentationScheduler, itk::LightObject);

  /** Get/Set methods */
  void SetMaxNumberOfIterations(const uint32_t niter){ m_MaxNumberOfIterations = niter;}
  void SetStartingNumberOfIterations(const uint32_t nbiter){ m_StartingNumberOfIterations = nbiter; }
  void SetPartialNumberOfIterations(const uint32_t nbiter){ m_PartialNumberOfIterations = nbiter; }
  void SetThreshold(const float thresh){ m_Threshold = thresh; }
  void SetSpectralWeight(const float sw){ m_SpectralWeight = sw; }
  void SetShapeWeight(const float spaw){ m_ShapeWeight = spaw; }
  void SetBandWeights(const std::vector<float> bw){ m_BandWeights = bw; }
  void SetAggregateGraphs(bool AggregateGraphs){m_AggregateGraphs = AggregateGraphs;};

protected:

  /** Constructor */
  LSBaatzSegmentationScheduler();

  /** Destructor */
  ~LSBaatzSegmentationScheduler();

  /**Create Filter*/
  itk::SmartPointer<BaatzSegmentationFilterType> CreateFilter();

  /** Generation method */
  virtual void GenerateData();

  virtual void ComputePaddingValue();

  virtual void NoTilingExecution();

  virtual void TilingExecution();

private:

	enum SegState{ 
		PARTIAL_SEG = 0, // Next step is a partial segmentation
		SEG_OVER_NO_AGGR, // Next step is to do nothing, segmentation is over and graphs cannot be merged
		AGGR_AND_SEG, // Next step is merging the graphs and achieve the segmentation
		NO_AGGR_AND_SEG,//No aggregation but keep segmenting
		AGGR_NO_SEG, // Next step is merging the graphs. 
		SEG_STATE_UNDEF
	};

  void RescaleGraph(ProcessingTile& tile);
  void AchieveSegmentation();
  void FinalGraphAgregation();
  void PartialSegmentation();
  void PartialSegmentation(uint32_t numberIterations);
  void ExtractStabilityMargins();
  void AggregateStabilityMargins();
  void RunPartialSegmentation(unsigned long int& accumulatedMemory, int& fusionSum);
  void NoAggregationAndPartialSegmentation();

  SegState FirstPartialSegmentation();

  /** Reachable member attributes */

  // Global maximum number of iterations to apply on the image.
  uint32_t m_MaxNumberOfIterations; 

  // Number of iterations for the first partial segmentation
  uint32_t m_StartingNumberOfIterations;

  // Number of iterations of the next partial segmentations
  uint32_t m_PartialNumberOfIterations;

  // Current number of iterations executed
  uint32_t m_CurrentNumberOfIterations;

  // Threshold to determine if two adjacent nodes are similar
  float m_Threshold;

  // Relative importance given to the spectral information
  float m_SpectralWeight;

  // Relative importance given to geometrical information
  float m_ShapeWeight;

  // Relative importance of each spectral band
  std::vector< float > m_BandWeights;

  /** Internal member attributes */

  // The serialized stability margin extracted from the current graph
  std::vector< char > m_SerializedStabilityMargin;

  // The maximum number of bytes to send via MPI
  unsigned long int m_MaxNumberOfBytes;

  // Flag to activate or deactivate the aggreagation of graph
  bool m_AggregateGraphs;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaLSBaatzSegmentationScheduler.txx"
#endif
