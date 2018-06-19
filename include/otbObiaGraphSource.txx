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
#ifndef otbObiaGraphSource_txx
#define otbObiaGraphSource_txx
#include "otbObiaGraphSource.h"

namespace otb
{
namespace obia
{

template< typename TOutputGraph >
GraphSource< TOutputGraph >
::GraphSource()
{
    // Create the output. We use static_cast<> here because we know the default
    // output must be of type TOutputGraph
    typename TOutputGraph::Pointer output =
    static_cast< TOutputGraph * >( this->MakeOutput(0).GetPointer() );
    this->ProcessObject::SetNumberOfRequiredOutputs(1);
    this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

    // Set the default behavior of an image source to NOT release its
    // output bulk data prior to GenerateData() in case that bulk data
    // can be reused (an thus avoid a costly deallocate/allocate cycle).
    this->ReleaseDataBeforeUpdateFlagOff();
}

/**
 *
 */
template< typename TOutputGraph >
itk::ProcessObject::DataObjectPointer
GraphSource< TOutputGraph >
::MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType)
{
    return TOutputGraph::New().GetPointer();
}


/**
 *
 */
template< typename TOutputGraph >
itk::ProcessObject::DataObjectPointer
GraphSource< TOutputGraph >
::MakeOutput(const itk::ProcessObject::DataObjectIdentifierType &)
{
    return TOutputGraph::New().GetPointer();
}

template< typename TOutputGraph >
typename GraphSource< TOutputGraph >::OutputGraphType*
GraphSource< TOutputGraph >
::GetOutput()
{
    // we assume that the first output is of the templated type
    return itkDynamicCastInDebugMode< OutputGraphType * >( this->GetPrimaryOutput() );
}

template< typename TOutputGraph >
const typename GraphSource< TOutputGraph >::OutputGraphType*
GraphSource< TOutputGraph >
::GetOutput() const
{
    // we assume that the first output is of the templated type
    return itkDynamicCastInDebugMode< const OutputGraphType * >( this->GetPrimaryOutput() );
}

/**
 *
 */
template< typename TOutputGraph >
typename GraphSource< TOutputGraph >::OutputGraphType*
GraphSource< TOutputGraph >
::GetOutput(unsigned int idx)
{
    TOutputGraph *out = dynamic_cast< TOutputGraph * >
                        ( this->ProcessObject::GetOutput(idx) );

    if ( out == ITK_NULLPTR && this->ProcessObject::GetOutput(idx) != ITK_NULLPTR )
    {
        itkWarningMacro (<< "Unable to convert output number " << idx << " to type " <<  typeid( OutputGraphType ).name () );
    }
    return out;
}


/**
 *
 */
template< typename TOutputGraph >
void
GraphSource< TOutputGraph >
::GraftOutput(itk::DataObject *graft)
{
    this->GraftNthOutput(0, graft);
}

/**
 *
 */
template< typename TOutputGraph >
void
GraphSource< TOutputGraph >
::GraftOutput(const DataObjectIdentifierType & key, itk::DataObject *graft)
{
    if ( !graft )
    {
        itkExceptionMacro(<< "Requested to graft output that is a ITK_NULLPTR pointer");
    }

    // we use the process object method since all out output may not be
    // of the same type
    itk::DataObject *output = this->ProcessObject::GetOutput(key);

    // Call GraftImage to copy meta-information, regions, and the pixel container
    output->Graft(graft);
}

/**
 *
 */
template< typename TOutputGraph >
void
GraphSource< TOutputGraph >
::GraftNthOutput(unsigned int idx, itk::DataObject *graft)
{
    if ( idx >= this->GetNumberOfIndexedOutputs() )
    {
        itkExceptionMacro(<< "Requested to graft output " << idx
                          << " but this filter only has " << this->GetNumberOfIndexedOutputs() << " indexed Outputs.");
    }
    this->GraftOutput( this->MakeNameFromOutputIndex(idx), graft );
}


} // end of namespace obia
} // end of namespace otb

#endif
