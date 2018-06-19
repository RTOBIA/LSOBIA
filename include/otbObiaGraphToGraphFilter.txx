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
#ifndef otbObiaGraphToGraphFilter_txx
#define otbObiaGraphToGraphFilter_txx
#include "otbObiaGraphToGraphFilter.h"

#include "otbMPIConfig.h"

namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputGraph >
GraphToGraphFilter<TInputGraph, TOutputGraph>
::GraphToGraphFilter()
{
    // Modify superclass default values, can be overridden by subclasses
    this->SetNumberOfRequiredInputs(1);
}

template< typename TInputGraph, typename TOutputGraph >
GraphToGraphFilter<TInputGraph, TOutputGraph>
::~GraphToGraphFilter()
{
}

template< typename TInputGraph, typename TOutputGraph >
void 
GraphToGraphFilter<TInputGraph, TOutputGraph>
::SetInput(const InputGraphType *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( 0,
                                    const_cast< InputGraphType * >( input ) );
}

template< typename TInputGraph, typename TOutputGraph >
void 
GraphToGraphFilter<TInputGraph, TOutputGraph>
::SetInput(unsigned int index, const InputGraphType *input)
{
    // Process object is not const-correct so the const_cast is required here
    this->itk::ProcessObject::SetNthInput( index,
                                    const_cast< InputGraphType * >( input ) );
}

/**
 *
 */
template< typename TInputGraph, typename TOutputGraph >
const typename GraphToGraphFilter<TInputGraph, TOutputGraph>::InputGraphType *
GraphToGraphFilter<TInputGraph, TOutputGraph>
::GetInput() const
{
    return itkDynamicCastInDebugMode< const InputGraphType * >( this->GetPrimaryInput() );
}

/**
 *
 */
template< typename TInputGraph, typename TOutputGraph >
const typename GraphToGraphFilter<TInputGraph, TOutputGraph>::InputGraphType *
GraphToGraphFilter<TInputGraph, TOutputGraph>
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

template< typename TInputGraph, typename TOutputGraph >
typename GraphToGraphFilter<TInputGraph, TOutputGraph>::OutputGraphType *
GraphToGraphFilter<TInputGraph, TOutputGraph>
::GetOutputByMove()
{
    InputGraphType * inputGraphPtr = dynamic_cast< InputGraphType * >( this->itk::ProcessObject::GetInput(0));
    auto outputGraph = this->GetOutput();
    outputGraph->GraftGraphByMove(inputGraphPtr);

    //Reset
    inputGraphPtr->Reset();

    return outputGraph;
}

template< typename TInputGraph, typename TOutputGraph >
void
GraphToGraphFilter<TInputGraph, TOutputGraph>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // end of namespace obia
} // end of namespace otb

#endif
