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
