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
#ifndef otbObiaPolygonizeFilter_txx
#define otbObiaPolygonizeFilter_txx

namespace otb
{
namespace obia
{

template <class TGraphType>
PolygonizeFilter<TGraphType>
::PolygonizeFilter()
{
	m_GraphToVectorFilter = GraphToVectorFilterType::New();
	m_SimplifyVectorFilter = SimplifyVectorFilterType::New(); 
	
	m_GraphToVectorFilter->SetXshift(0);
	m_GraphToVectorFilter->SetYshift(0);
	
	auto simplifyFunc = new otb::obia::DouglasPeukerFunc();
	simplifyFunc->SetTolerance(1.0);
	
	
	m_SimplifyVectorFilter->SetSimplifyFunc(simplifyFunc);
	m_SimplifyVectorFilter->SetLayerName(cleanedLayerName);
}

template <class TGraphType>
void
PolygonizeFilter<TGraphType>::
GenerateData()
{
	m_GraphToVectorFilter->SetInput(this->GetInput());
	// TODO: Check why we have to update for each filter?
	m_GraphToVectorFilter->Update();
	
	m_SimplifyVectorFilter->SetInput(m_GraphToVectorFilter->GetOutput());
	
	m_SimplifyVectorFilter->GraftOutput(this->GetOutput());
	m_SimplifyVectorFilter->Update();
	this->GraftOutput(m_SimplifyVectorFilter->GetOutput());
}

} /* namespace obia */
} /* namespace otb */

#endif
