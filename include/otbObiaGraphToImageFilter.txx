/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

    const double origin[2] = { graph->GetOriginX(),  graph->GetOriginY()};
    image->SetOrigin(origin);
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
