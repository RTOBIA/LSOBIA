#ifndef otbObiaImageToGraphFilter_txx
#define otbObiaImageToGraphFilter_txx
#include "otbObiaImageToGraphFilter.h"

namespace otb
{
namespace obia
{

template< typename TInputImage, typename TOutputGraph >
ImageToGraphFilter<TInputImage, TOutputGraph>::
ImageToGraphFilter():
m_ProcessNoData(false),
m_NoDataValue(0)
{
      // Modify superclass default values, can be overridden by subclasses
      this->SetNumberOfRequiredInputs(1);
}

template< typename TInputImage, typename TOutputGraph >
ImageToGraphFilter<TInputImage, TOutputGraph>::
~ImageToGraphFilter()
{
}

template< typename TInputImage, typename TOutputGraph >
void 
ImageToGraphFilter<TInputImage, TOutputGraph>::
SetInput(const InputImageType *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( 0,
                                      const_cast< InputImageType * >( input ) );
}

template< typename TInputImage, typename TOutputGraph >
void 
ImageToGraphFilter<TInputImage, TOutputGraph>::
SetInput(unsigned int index, const TInputImage *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( index,
                                      const_cast< TInputImage * >( input ) );
}

/**
 *
 */
template< typename TInputImage, typename TOutputGraph >
const typename ImageToGraphFilter<TInputImage, TOutputGraph>::InputImageType *
ImageToGraphFilter<TInputImage, TOutputGraph>
::GetInput() const
{
    return itkDynamicCastInDebugMode< const TInputImage * >( this->GetPrimaryInput() );
}

/**
 *
 */
template< typename TInputImage, typename TOutputGraph >
const typename ImageToGraphFilter<TInputImage, TOutputGraph>::InputImageType *
ImageToGraphFilter<TInputImage, TOutputGraph>
::GetInput(unsigned int idx) const
{
    const TInputImage *in = dynamic_cast< const TInputImage * >
      ( this->itk::ProcessObject::GetInput(idx) );

    if ( in == ITK_NULLPTR && this->itk::ProcessObject::GetInput(idx) != ITK_NULLPTR )
    {
        itkWarningMacro (<< "Unable to convert input number " << idx << " to type " <<  typeid( InputImageType ).name () );
    }
    return in;
}

template< typename TInputImage, typename TOutputGraph >
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // end of namespace obia
} // end of namespace otb

#endif
