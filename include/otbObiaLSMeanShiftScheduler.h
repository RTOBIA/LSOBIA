#ifndef otbObiaLSMeanShiftScheduler_h
#define otbObiaLSMeanShiftScheduler_h
#include "otbObiaLSImageToGraphScheduler.h"
#include "otbObiaLabelImageToGraphFilter.h"
#include "otbImage.h"
#include "otbMeanShiftSmoothingImageFilter.h"
#include "otbConnectedComponentMuParserFunctor.h"
#include "itkConnectedComponentFunctorImageFilter.h"

namespace otb
{
namespace obia
{

template< typename TInputImage, typename TLabelPixel>
class LSMeanShiftScheduler : public LSImageToGraphScheduler< TInputImage,
    typename LabelImageToGraphFilter<TLabelPixel>::OutputGraphType >
{
    public:

    /** Some convenient class alias */
    using InputImageType = TInputImage;
    using InputInternalPixelType = typename InputImageType::InternalPixelType;
    using InputPixelType = typename InputImageType::PixelType;
    using InputImageFileReaderType = ImageFileReader<InputImageType>;
    using MultiChannelExtractROIFilterType = MultiChannelExtractROI<InputInternalPixelType, InputInternalPixelType>;
    using MeanShiftSmoothingImageFilterType = MeanShiftSmoothingImageFilter<InputImageType, InputImageType>;
    using LabelPixelType = TLabelPixel;
    using LabelImageType = Image<LabelPixelType>;
    using CCFunctorType = otb::Functor::ConnectedComponentMuParserFunctor<InputPixelType>;
    using CCFilterType = itk::ConnectedComponentFunctorImageFilter<InputImageType, LabelImageType, CCFunctorType, otb::Image<unsigned int> >;
    using LabelImageToGraphFilterType = LabelImageToGraphFilter<LabelPixelType>;
    using OutputGraphType = typename LabelImageToGraphFilterType::OutputGraphType;
    using OutputGraphPointerType = typename OutputGraphType::Pointer;
    using NodeType = typename OutputGraphType::NodeType;
    using GraphOperationsType = GraphOperations<OutputGraphType>;

    /** Standard class alias */
    using Self            = LSMeanShiftScheduler;    
    using SuperClass   = LSImageToGraphScheduler<InputImageType, OutputGraphType>;
    using Pointer      = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >; 
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(LSMeanShiftScheduler, LSImageToGraphScheduler);

    inline void SetMaxNumberOfIterations(const unsigned int niter){ m_MaxNumberOfIterations = niter; }
    inline void SetSpatialBandWidth(const unsigned int spatialr){ m_SpatialBandWidth = spatialr; }
    inline void SetSpectralRangeBandWidth(const float spectralr){ m_SpectralRangeBandWidth = spectralr; }
    inline void SetThreshold(const float thresh){ m_Threshold = thresh; }
    inline void SetSpectralRangeRamp(const float ranger){ m_SpectralRangeRamp = ranger; }


protected:
    
    LSMeanShiftScheduler();
    virtual ~LSMeanShiftScheduler();

    /** Generation method */
    virtual void GenerateData();

    /** 
    Compute the stability margin of the Mean-Shift algorithm using
    the Connected Component algorithm. (Refer to the publication of
    David Youssefi and Julien Michel).
     */
    virtual void ComputePaddingValue();

    /** Mean Shift -> Connected Component -> Label To Graph */
    void Segment();

    /** Achieve the construction of the nodes located at the borders of the graph */
    void PostProcessing();

private:

    void ComputeExtractionParametersForCC(ProcessingTile& tile,
                                          uint32_t& startX,
                                          uint32_t& startY);

    void RescaleGraph(ProcessingTile& tile);

    void SerializeListOfBorderNodes(const ProcessingTile& tile);

    void FillSharedBuffer();

    void AchieveBorderNodeConstruction();

    /**\brief Aggregate all graphs to the master graph (rank 0)
     * */
    void FinalGraphAgregation();

    /** Algorithm iterative scheme will stop if convergence hasn't been reached after the maximum number of iterations. */
    unsigned int m_MaxNumberOfIterations;

    /** Spatial radius of the neighborhood. */
    unsigned int m_SpatialBandWidth;

    /** Range radius defining the radius (expressed in radiometry unit) in the multi-spectral space. */
    float m_SpectralRangeBandWidth;

    /** Algorithm iterative scheme will stop if mean-shift vector is below this threshold or if iteration number reached maximum number of iterations. */
    float m_Threshold;

    /** This coefficient makes dependent the ranger of the colorimetry of the filtered pixel : y = rangeramp*x+ranger. */
    float m_SpectralRangeRamp;

    // The serialized list of border nodes.
    std::vector< char > m_SerializedListOfBorderNodes;

    // The maximum number of bytes to send via MPI
    unsigned long int m_MaxNumberOfBytes;

    // The shared buffer that can be accessed by other processors
    std::vector< char > m_SharedBuffer;

}; 

} // end of namespace obia
} // end of namespace otb
#include "otbObiaLSMeanShiftScheduler.txx"
#endif
