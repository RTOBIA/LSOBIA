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
#ifndef otbObiaGraphSource_h
#define otbObiaGraphSource_h
#include "itkProcessObject.h"
#include "otbObiaGraph.h"

namespace otb
{
namespace obia
{

/** \class GraphSource
 *    \brief Base class for all process objects that output graph data.
 *
 * GraphSource is the base class for all process objects that output
 * graph data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the adjacent graph.
 */
template< typename TOutputGraph >
class GraphSource : public itk::ProcessObject
{

public:

    /** Standard class alias. */
    using Self         = GraphSource;
    using Superclass   = itk::ProcessObject;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self >;

    /** Smart Pointer type to a DataObject. */
    using DataObjectPointer = itk::DataObject::Pointer;

    using DataObjectIdentifierType       = Superclass::DataObjectIdentifierType;
    using DataObjectPointerArraySizeType = Superclass::DataObjectPointerArraySizeType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(GraphSource, ProcessObject);

    /** Some convenient alias */
    using OutputGraphType = TOutputGraph;

    OutputGraphType * GetOutput();
    const OutputGraphType * GetOutput() const;
    OutputGraphType * GetOutput(unsigned int idx);
    OutputGraphType* GetOutputByMove();

    virtual itk::ProcessObject::DataObjectPointer 
          MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

    virtual itk::ProcessObject::DataObjectPointer 
          MakeOutput(const itk::ProcessObject::DataObjectIdentifierType &) ITK_OVERRIDE;

    virtual void GraftOutput(itk::DataObject *output);
    virtual void GraftOutput(const DataObjectIdentifierType & key, itk::DataObject *output);
    virtual void GraftNthOutput(unsigned int idx, itk::DataObject *output);

protected:

    GraphSource();
    virtual ~GraphSource() {}

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphSource.txx"
#endif
