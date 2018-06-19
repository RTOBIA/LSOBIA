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
#ifndef otbObiaComputeAttributesFilter_txx
#define otbObiaComputeAttributesFilter_txx

#include <otbObiaComputeAttributesFilter.h>
#include "itkDefaultConvertPixelTraits.h"
#include "otbObiaVectorOperations.h"

namespace otb
{
namespace obia
{

/** Constructor */
template <class TInputImage, class TMaskImage>
ComputeAttributesFilter<TInputImage,TMaskImage>
::ComputeAttributesFilter() :   m_NumberOfBands(1),
								m_InputImage(nullptr)
{
    // Modify superclass default values, can be overridden by subclasses
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfThreads(1);

}

/**
 * Reset the persistent data of the filter.
 */
template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::Reset()
{

}
/**
 * Synthesize the persistent data of the filter.
 */
template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::Synthetize()
{

}
/** Generate data should thread over */
template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::GenerateData()
{
	//Initialize output layer
	this->InitializeOutput();

	//Allocate output memory
	//this->AllocateOutputs();

	//this->BeforeThreadedGenerateData();

	// Split the data into in-memory layers
	//this->DispatchInputVectors();

	//Compute all attributs
	this->ComputeAllAttributes();

	//Write feature
	std::cout << "------------------- WRITING OUTPUT FILE -----------------" << std::endl;
	const OGRDataSourceType* ogrDs = this->GetOGRData();
//	std::stringstream ss;
//	ss << this->m_OutputDir << m_OutputFilename << ".gml";
	VectorOperations::WriteOGRDataSource(m_OutputDs, m_OutputFilename, "");
}

template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::InitializeOutput()
{
	std::cout << "-------------- INITIALIZE OUTPUT -------------" << std::endl;
	//Create a new layer
	const OGRDataSourceType* ogrDs = this->GetOGRData();

	//initialize output ds pointer
	m_OutputDs = OGRDataSourceType::New();

	//Create all new fields
	for(unsigned int k = 0; k < m_Attributes.size(); ++k)
	{
		GenericAttributeType* currentAttribut = m_Attributes[k];
		this->CreateAdditionalField(currentAttribut->GetFieldName(), currentAttribut->GetFieldType());
	}

	//Call method for adding fields
	this->InitializeOutputDataSource(const_cast<OGRDataSourceType* >(ogrDs), m_OutputDs);

	//Set primary output
	this->SetPrimaryOutput(m_OutputDs);
}

template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::ComputeAllAttributes()
{
	std::cout << "-------- COMPUTING ALL ATTRIBUTS -------------------------" << std::endl;

	// Retrieve inputs
	m_InputImage = const_cast<InputImageType*>(this->GetInput());
	const OGRDataSourceType* ogrDs = this->GetOGRData();

	//Set origin to 0
	typename InputImageType::PointType newOrigin;
	newOrigin.Fill(0.0);
	m_InputImage->SetOrigin( newOrigin );
//
//	double origin[2] = {0.0, 0.0};
//	m_InputImage->SetOrigin(origin);

	double spacing[2] = {1, 1};
	m_InputImage->SetSpacing(spacing);

	m_NumberOfBands = m_InputImage->GetNumberOfComponentsPerPixel();
	std::cout << "Number of bands = " << m_NumberOfBands << std::endl;

	RegionType requestedRegion = m_InputImage->GetLargestPossibleRegion();


	itk::Indent indent;
	// Loop across the features in the layer (filtered by requested region in BeforeTGD already)
	OGRLayerType layer = ogrDs->GetLayer(m_InputLayerName);
	ogr::Layer::const_iterator featIt = layer.begin();


	//TODO : Do we need to use dynamic thread id?
	itk::ThreadIdType threadId = 0;

	for(; featIt!= layer.end(); ++featIt)
	{
		//Current feature
		OGRFeatureType currentFeature = *featIt;
		//std::cout << "Feature " << currentFeature.GetFID() << std::endl;

		//Copy the feature to the outputlayer
		OGRLayer* outputLayer = m_OutputDs->ogr().GetLayerByName(this->GetOutLayerName().c_str());

		//Generate feature to add fields
		OGRFeatureType updatedFeature(*(outputLayer->GetLayerDefn()));
		updatedFeature.SetFrom(currentFeature, true);

		//Loop attributs
		for(unsigned int k = 0; k < m_Attributes.size(); k++)
		{
			//Clear m_Samples
			m_Samples.clear();

			//Get current attribut class
			GenericAttributeType* currentAttribut = m_Attributes[k];

			//Use extract shape from the attribut method
			OGRFeatureType requiredFeature = currentAttribut->GenerateFeature(currentFeature);

			//Use this shape for compute region and extract pixel
			RegionType consideredRegion = this->FeatureBoundingRegion(m_InputImage, featIt);
			bool regionNotEmpty = consideredRegion.Crop(requestedRegion);

			if(regionNotEmpty)
			{
				//std::cout << "Region not empty : " << consideredRegion << std::endl;
				//std::cout << "Number of pixels " << consideredRegion.GetNumberOfPixels() << std::endl;
				this->ExploreGeometry(currentFeature, currentFeature.ogr().GetGeometryRef(),consideredRegion, threadId);

				//Give these pixel to the attribut method for each band
				double tmp = 0.0;
				for(unsigned int b = 0; b < m_NumberOfBands; b++)
				{
					tmp += currentAttribut->ComputeAttribut(m_Samples[b]);
				}
				//Compute the mean over all the bands
				tmp /= m_NumberOfBands;

				//Update feature
				updatedFeature.ogr().SetField(currentAttribut->GetFieldName().c_str(), tmp);
			}
		}

		//Update layer
		OGRErr errCreate = outputLayer->CreateFeature(&(updatedFeature.ogr()));
		if(errCreate != 0)
		{
			std::cerr << "Error when creating feature " << updatedFeature.ogr().GetFID() << " in " << __func__ << std::endl;
		}
	}

}


/** Generic method called for each matching pixel position*/
template <class TInputImage, class TMaskImage>
void
ComputeAttributesFilter<TInputImage,TMaskImage>
::ProcessSample(const ogr::Feature& feature,
			    typename InputImageType::IndexType& imgIndex,
			    typename InputImageType::PointType& imgPoint,
			    itk::ThreadIdType& threadid)
{
//	std::cout << "----------------------------------------------" << std::endl;
//	std::cout << "Process sample : " << imgPoint << std::endl;
//	std::cout << "Image index    : " << imgIndex << std::endl;
//	std::cout << "----------------------------------------------" << std::endl;
	//m_Samples.push_back()
	typename InputImageType::PixelType pixel = m_InputImage->GetPixel(imgIndex);

	for (unsigned int b = 0 ; b < m_NumberOfBands ; ++b)
	{
		//std::cout << "Pixel " << i << " = "<< pixel[b] << std::endl;
		m_Samples[b].push_back(pixel[b]);
	}
}


}//end obia

}//end otb

#endif
