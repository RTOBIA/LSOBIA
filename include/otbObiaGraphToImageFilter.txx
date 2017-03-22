#ifndef otbObiaGraphToImageFilter_txx
#define otbObiaGraphToImageFilter_txx
#include "otbObiaGraphToImageFilter.h"

namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputImage >
GraphToImageFilter<TInputGraph, TOutputImage>::
GraphToImageFilter()
{
    // Modify superclass default values, can be overridden by subclasses
    this->SetNumberOfRequiredInputs(1);
}

template< typename TInputGraph, typename TOutputImage >
GraphToImageFilter<TInputGraph, TOutputImage>::
~GraphToImageFilter()
{
}

template< typename TInputGraph, typename TOutputImage >
void 
GraphToImageFilter<TInputGraph, TOutputImage>::
SetInput(const InputGraphType *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( 0,
                                    const_cast< InputGraphType * >( input ) );
}

template< typename TInputGraph, typename TOutputImage>
void
GraphToImageFilter<TInputGraph, TOutputImage>::
UpdateOutputInformation()
{
    Superclass::UpdateOutputInformation();
    auto graph = this->GetInput();

    typename OutputImageType::IndexType index;
    typename OutputImageType::SizeType size;
    typename OutputImageType::RegionType region;

    index[0] = 0; index[1] = 0;
    size[0] = graph->GetImageWidth();
    size[1] = graph->GetImageHeight();

    region.SetIndex(index);
    region.SetSize(size);

    auto image = this->GetOutput();
    image->SetRegions(region);
    image->SetProjectionRef(graph->GetProjectionRef());
    image->SetNumberOfComponentsPerPixel(graph->GetNumberOfSpectralBands());
}


template< typename TInputGraph, typename TOutputImage >
void GraphToImageFilter<TInputGraph, TOutputImage>::GenerateOutputInformation()
{
}


template< typename TInputGraph, typename TOutputImage >
void 
GraphToImageFilter<TInputGraph, TOutputImage>::
SetInput(unsigned int index, const InputGraphType *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( index,
                                    const_cast< InputGraphType * >( input ) );
}

/**
 *
 */
template< typename TInputGraph, typename TOutputImage >
const typename GraphToImageFilter<TInputGraph, TOutputImage>::InputGraphType *
GraphToImageFilter<TInputGraph, TOutputImage>
::GetInput() const
{
    return itkDynamicCastInDebugMode< const InputGraphType * >( this->GetPrimaryInput() );
}

/**
 *
 */
template< typename TInputGraph, typename TOutputImage >
const typename GraphToImageFilter<TInputGraph, TOutputImage>::InputGraphType *
GraphToImageFilter<TInputGraph, TOutputImage>
::GetInput(unsigned int idx) const
{
  const InputGraphType *in = dynamic_cast< const InputGraphType * >
    ( this->itk::ProcessObject::GetInput(idx) );

    if ( in == ITK_NULLPTR && this->itk::ProcessObject::GetInput(idx) != ITK_NULLPTR )
    {
        itkWarningMacro (<< "Unable to convert input number " << idx << " to type " <<  typeid( InputGraphType ).name () );
    }
    return in;
}

template< typename TInputGraph, typename TOutputImage >
void
GraphToImageFilter<TInputGraph, TOutputImage>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // end of namespace obia
} // end of namespace otb
#endif
