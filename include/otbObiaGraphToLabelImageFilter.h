#ifndef otbObiaGraphToLabelImageFilter_h
#define otbObiaGraphToLabelImageFilter_h
#include "itkImageRegionIterator.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "otbObiaGraphToImageFilter.h"


namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputImage>
class GraphToLabelImageFilter : public GraphToImageFilter<TInputGraph, TOutputImage>
{
public:

    /** Some convenient alias */
    using Self                       = GraphToLabelImageFilter;
    using Superclass                 = itk::ImageSource<TOutputImage>;
    using Pointer                    = itk::SmartPointer<Self>;
    using ConstPointer               = itk::SmartPointer< const Self >;
    using InputGraphType             = TInputGraph;
    using InputGraphPointerType  	 = typename TInputGraph::Pointer;
    using NodeType                	 = typename InputGraphType::NodeType;
    using OutputImageType         	 = TOutputImage;
    using OutputImagePointerType  	 = typename TOutputImage::Pointer;
    using OutputImageIteratorType  	 = itk::ImageRegionIterator<OutputImageType>;
    using FillholeFilterType      	 = itk::GrayscaleFillholeImageFilter<OutputImageType,OutputImageType>;

    itkNewMacro(Self);

    //Get LUT
    std::vector<typename OutputImageType::InternalPixelType> GetReverseLut(){return m_ReverseLut;};

protected:

    virtual void GenerateData();

    //Keep LUT to get correspondance between label and node id
    std::vector<typename OutputImageType::InternalPixelType> m_ReverseLut;
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraphToLabelImageFilter.txx"
#endif
