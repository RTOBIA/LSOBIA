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
#ifndef otbObiaMeanAttribute_txx
#define otbObiaMeanAttribute_txx

#include <otbObiaMeanAttribute.h>

namespace otb
{
namespace obia
{

template< class TInputImage >
MeanAttribute<TInputImage>
::MeanAttribute()
{
	this->m_FieldName 	 = "Mean";
	this->m_AttributName = "Mean Attribute";
	this->m_FieldType	 = OFTReal;
}

template< class TInputImage >
MeanAttribute<TInputImage>
::~MeanAttribute()
{

}

/**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
template< class TInputImage >
typename MeanAttribute<TInputImage>::OGRFeatureType
MeanAttribute<TInputImage>
::GenerateFeature(const OGRFeatureType inputFeature)
{
	//Create a new feature (clone of input feature because for the mean, we need the full geometry)
	OGRFeatureType requiredFeature(inputFeature);

	return requiredFeature;
}

/**Compute the attribut*/
template< class TInputImage >
double
MeanAttribute<TInputImage>
::ComputeAttribut(std::vector<PixelType> samples)
{
	double meanValue = 0.0;
	//Loop each samples
	for(unsigned int k = 0; k < samples.size(); k++)
	{
		PixelType sample = samples[k];
		meanValue += static_cast<double>(sample);
	}

	if(samples.size() != 0)
	{
		meanValue = meanValue/samples.size();
	}

	return meanValue;
}

}//end obia
}//end otb

#endif
