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
#include "otbObiaComputeAttributesFilter.h"
#include "otbObiaComputeAttributesFilter.h"
#include "otbObiaGenericAttribute.h"
#include "otbObiaMeanAttribute.h"
#include "otbObiaBorderAttribute.h"
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbObiaConstExpr.h"
#include "otbImageFileReader.h"

#define NUM_ELEMENT 4

using InputImageType = otb::VectorImage<float, 2>;
using PixelType      = InputImageType::PixelType;
std::vector<otb::obia::GenericAttribute<InputImageType>*> CreateAttributes();

bool isOn(const std::string & str)
{
	if (str == std::string("on"))
	{
		return true;
	}
	return false;
}

int otbLSOBIAComputeAttributesTest(int argc, char *argv[])
{
    if(argc < 3)
    {

        std::cerr << "Usage " << argv[0] << ": \n"
            << "Argument 1: " << "[input vector file]\n"
			<< "Argument 2: " << "[input image file]\n"
            << "Argument 3: " << "[temporary directory to store intermediate files for a node]\n"
            << "Argument 4: " << "[output directory path (must me already created)]\n"
			<< "Argument 5: " << "[Output file name]\n"
			<< std::endl;
        return 1;
    }


    /* Global parameters */
	std::string vectorFile	= argv[1];
	std::string imagePath 	= argv[2];
	std::string tmpDir 		= argv[3];
	std::string outDir 		= argv[4];
	std::string gmlFile 	= argv[5];

	//Read image
	//Image reader
	auto imgReader = otb::ImageFileReader<InputImageType>::New();
	imgReader->SetFileName(imagePath);
	imgReader->Update();

	//Compute attributes
	using OGRDataSourceType = otb::ogr::DataSource;
	using ComputeAttributesFilterType = otb::obia::ComputeAttributesFilter<InputImageType>;
	auto computeAttributes = ComputeAttributesFilterType::New();

	OGRDataSourceType::Pointer outputDs = otb::ogr::DataSource::New(vectorFile, OGRDataSourceType::Modes::Read); //default 2nd argument is read
	std::cout << "Number of layers = " << outputDs->GetLayersCount() << std::endl;
	std::cout << "Layer = " << outputDs->GetLayer(0).GetName() << std::endl;
	computeAttributes->SetOGRData(outputDs);
	computeAttributes->SetInput(imgReader->GetOutput());
	computeAttributes->SetInputLayerName(otb::obia::reconstructedLayerName);
	computeAttributes->SetLayerIndex(0);
	computeAttributes->SetOutLayerName(otb::obia::attributesLayerName);
	computeAttributes->SetFieldName(otb::obia::startingCoordsFieldName);
	computeAttributes->SetAttributes(CreateAttributes());
	computeAttributes->SetOutputDir(outDir);
	computeAttributes->SetOutputFilename(gmlFile);
	computeAttributes->Update();

    std::cout << "SUCCESS" << std::endl;
    return EXIT_SUCCESS;
}

std::vector<otb::obia::GenericAttribute<InputImageType>*> CreateAttributes()
{
	using GenericAttributeType 	= otb::obia::GenericAttribute<InputImageType>;
	using MeanAttributeType 	= otb::obia::MeanAttribute<InputImageType>;
	using BorderAttributeType 	= otb::obia::BorderAttribute<InputImageType>;

	//Vector of attributes
	std::vector<otb::obia::GenericAttribute<InputImageType>*> attributes;

	//Create first attribute
	GenericAttributeType* meanAttribute = new MeanAttributeType();
	attributes.push_back(meanAttribute);

	//Create second attribute
	GenericAttributeType* borderAttribute = new BorderAttributeType(2.0);
	attributes.push_back(borderAttribute);


	return attributes;
}
