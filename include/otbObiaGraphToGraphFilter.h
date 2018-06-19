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
#ifndef otbObiaGraphToGraphFilter_h
#define otbObiaGraphToGraphFilter_h
#include "otbObiaGraphSource.h"

namespace otb
{
namespace obia
{

/** \class GraphToGraphFilter
 *    \brief Base class for all process objects that requires one input graph
     and that provides one output graph data.
 *
 */
template< typename TInputGraph, typename TOutputGraph >
class GraphToGraphFilter : public GraphSource< TOutputGraph >
{

public:

    /** Standard class alias */
    using Self         = GraphToGraphFilter;
    using Superclass   = GraphSource<TOutputGraph>;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self >;

    /** Run-time type information (and related methods). */
    itkTypeMacro(GraphToGraphFilter, GraphSource);
    itkNewMacro(Self);

    /** Some convenient typedefs. */
    using InputGraphType         = TInputGraph;
    using InputGraphPointer      = typename InputGraphType::Pointer;
    using InputGraphConstPointer = typename InputGraphType::ConstPointer;

    using OutputGraphType         = TOutputGraph;
    using OutputGraphPointer      = typename OutputGraphType::Pointer;
    using OutputGraphConstPointer = typename OutputGraphType::ConstPointer;

    virtual void SetInput(const InputGraphType *input);
    virtual void SetInput(unsigned int, const InputGraphType *image);

    const InputGraphType * GetInput() const;
    const InputGraphType * GetInput(unsigned int idx) const;

    OutputGraphType* GetOutputByMove();


protected:

    GraphToGraphFilter();
    ~GraphToGraphFilter();

    virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

private:

    GraphToGraphFilter(const Self &) =delete;
    void operator=(const Self &) =delete;

};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraphToGraphFilter.txx"
#endif
