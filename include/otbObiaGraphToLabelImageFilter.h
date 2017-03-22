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
    using InputGraphType           = TInputGraph;
    using InputGraphPointerType   = typename TInputGraph::Pointer;
    using NodeType                = typename InputGraphType::NodeType;
    using OutputImageType           = TOutputImage;
    using OutputImagePointerType  = typename TOutputImage::Pointer;
    using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;
    using FillholeFilterType      = itk::GrayscaleFillholeImageFilter<OutputImageType,OutputImageType>;

    itkNewMacro(Self);

protected:

    virtual void GenerateData();
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraphToLabelImageFilter.txx"
#endif
