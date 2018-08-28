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
#ifndef otbObiaPolygonizeFilter_h
#define otbObiaPolygonizeFilter_h
#include "otbObiaGraphToVectorFilter.h"
#include "otbObiaSimplifyVectorFilter.h"
#include "otbObiaDouglasPeukerSimplify.h"

namespace otb
{
namespace obia
{

/**Class Converts OBIA graph to OGR graph */
template <class TGraphType>
  class PolygonizeFilter :
    public GraphToVectorFilter<TGraphType>
  {
  public:
	/** Standard PolygonizeFilter alias */
	typedef PolygonizeFilter                  Self;
	typedef GraphToVectorFilter<TGraphType>   Superclass;
	typedef SmartPointer<Self>                Pointer;
	typedef SmartPointer<const Self>          ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(PolygonizeFilter, GraphToVectorFilter);

    /** Template parameters typedefs */
    using InputGraphType = TGraphType;

    /** Constructor */
    PolygonizeFilter();
    ~PolygonizeFilter();

    /** Does the real work. */
    virtual void GenerateData();

  private:
    typedef otb::obia::GraphToVectorFilter<InputGraphType> GraphToVectorFilterType;
    typedef otb::obia::SimplifyVectorFilter<otb::obia::DouglasPeukerFunc> SimplifyVectorFilterType;

    typename GraphToVectorFilterType::Pointer m_GraphToVectorFilter;
    typename SimplifyVectorFilterType::Pointer m_SimplifyVectorFilter;

  };
} /* namespace obia */
} /* namespace otb */

#include "otbObiaPolygonizeFilter.txx"
#endif
