#ifndef otbObiaLSBaatzSegmentationScheduler_h
#define otbObiaLSBaatzSegmentationScheduler_h
#include "otbObiaBaatzSegmentationFilter.h"
#include "otbObiaLSImageToGraphScheduler.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaGenericRegionMergingFilter.h"


/**
\file otbObiaLSBaatzSegmentationScheduler.h
\brief This file define the class baatz graph scheduler used to handle multi-processing of baatz filtering
*/
namespace otb
{
namespace obia
{

/**\class LSBaatzSegmentationScheduler otbObiaLSBaatzSegmentationScheduler.h
 * \brief Class handling multi-processing of baatz segmentation, sharing stability margins, etc ...*/
template< class TInputImage >
  class LSBaatzSegmentationScheduler : public LSImageToGraphScheduler< TInputImage,
                        typename ImageToBaatzGraphFilter<TInputImage>::OutputGraphType >
{
public:

    /** Standard class alias */
    using Self          = LSBaatzSegmentationScheduler;
    using SuperClass   = itk::Object;
    using Pointer      = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >; 

    /** Some convenient alias */
    using InputImageType                   = TInputImage;
    using InputInternalPixelType           = typename InputImageType::InternalPixelType;
    using InputImageReaderType             = ImageFileReader<InputImageType>;
    using MultiChannelExtractROIFilterType = MultiChannelExtractROI<InputInternalPixelType,InputInternalPixelType>;
    using ImageToBaatzGraphFilterType        = ImageToBaatzGraphFilter<InputImageType>;
    using GraphType                        = typename ImageToBaatzGraphFilterType::OutputGraphType;
    using NodeType                         = typename GraphType::NodeType;
    using GraphPointerType                 = typename GraphType::Pointer;
    using GraphOperationsType              = GraphOperations<GraphType>;
    using BaatzSegmentationFilterType         = GenericRegionMergingFilter<GraphType, GraphType,
                                                                                                    BaatzMergingCost<float, GraphType> ,
                                                                      BaatzHeuristic<GraphType>,
                                                                      BaatzUpdateAttribute<GraphType> >;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(LSBaatzSegmentationScheduler, itk::LightObject);

    /**\brief Set the maximum number of iterations (to prevent infinite looping)
     * \param: Max number of iterations */
    void SetMaxNumberOfIterations(const uint32_t niter){ m_MaxNumberOfIterations = niter;}

    /**\brief Set the starting number of iterations (usually 0)
     * \param: Starting number of iterations */
    void SetStartingNumberOfIterations(const uint32_t nbiter){ m_StartingNumberOfIterations = nbiter; }


    /**\brief Set the partial number of iterations. Used in BaatzSegmentationFilter to compute nPartialIteration
     * before next step. Required to compute number of adjacency layers.
     * \param: Partial number of iterations */
    void SetPartialNumberOfIterations(const uint32_t nbiter){ m_PartialNumberOfIterations = nbiter; }

    /**\brief Set the threshold for baatz segmentation criterion
	 * \param: Threshold*/
    void SetThreshold(const float thresh){ m_Threshold = thresh; }

    /**\brief Set the minimum decreasing for accumulated memory during baatz segmentation
	 * \param: decreasing*/
	void SetDecreasing(const float decreasing){ m_Decreasing = decreasing; }

    /**\brief Set the spectral weight, used for baatz segmentation criterion
	 * \param: Spectral weight*/
    void SetSpectralWeight(const float sw){ m_SpectralWeight = sw; }

    /**\brief Set the shape weight, used for baatz segmentation criterion
 	 * \param: Shape weight*/
    void SetShapeWeight(const float spaw){ m_ShapeWeight = spaw; }

    /**\brief Set the band weight, used for baatz segmentation criterion
 	 * \param: band weight*/
    void SetBandWeights(const std::vector<float> bw){ m_BandWeights = bw; }

    /**\brief Set the boolean to indicate if graphs need to be merged (to save some memory for example, or to
     * get only 1 graph as output)
 	 * \param: Boolean*/
    void SetAggregateGraphs(bool AggregateGraphs){m_AggregateGraphs = AggregateGraphs;};

protected:

    /** Constructor */
    LSBaatzSegmentationScheduler();

    /** Destructor */
    ~LSBaatzSegmentationScheduler();

    /**Create Filter*/
    itk::SmartPointer<BaatzSegmentationFilterType> CreateFilter();

    /**\brief Implementation of generate data. It checks if the process is with tiles or not*/
    virtual void GenerateData();

    /**\brief Compute required padding for stability margin*/
    virtual void ComputePaddingValue();

    /**\brief Process in case of no tiles. it only load the input image and run baatz filter segmentation.
     * No need to handle stability margin*/
    virtual void NoTilingExecution();


    /**\brief Process in case of tiles. It first execute some iterations before deciding if graphs need to be
     * aggregate, or if segmentation keep going on each graph.It will loop over all tiles, extract all required stability margins,
     * aggregate these margins, and then run filter of each aggregate graphs, and so on. Synchronization is required for
     * following steps:
     * - Extract stability margins
     * - Aggregate stability margins*/
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
    /**\brief Modify coordinate of the graph associated to the tile to be in a absolute reference (and not tile reference)
     * \param: Tile associated to graph
     * */
    void RescaleGraph(ProcessingTile& tile);

    /**\brief Achieve segmentation if graphs are aggregated and current process is rank 0
     * */
    void AchieveSegmentation();

    /**\brief Aggregate all graphs to the master graph (rank 0)
     * */
    void FinalGraphAgregation();

    /**\brief Run m_PartialNumberOfIterations on current graph. It will call BaatzGraphFilter
     * with m_PartialNumberOfIterations iteration
     * */
    void PartialSegmentation();

    /**\brief Extract stability margin of the current graph
     * */
    void ExtractStabilityMargins();

    /**\brief Aggregate all stability margins to the current graph
     * */
    void AggregateStabilityMargins();


    /**\brief Run partial segmentation
     * */
    void RunPartialSegmentation(unsigned long int& accumulatedMemory, int& fusionSum);

    /**\brief In case of no aggregation and we still doing partial segmentation
     * */
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

    // Minimum decreasing wanted for accumulated memory
    float m_Decreasing;

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
