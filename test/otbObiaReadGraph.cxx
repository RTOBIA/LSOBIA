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
#include "otbObiaGraphFileReader.h"
#include "otbImage.h"
#include "otbObiaImageToBaatzGraphFilter.txx"
#include "otbObiaGraphOperations.txx"
#include "otbObiaGraphToLabelImageFilter.txx"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImageFileWriter.h"


int otbObiaReadGraph(int argc, char * argv[])
{
	if(argc < 2)
	{
	  std::cerr << "Usage " << argv[0] << ": \n"
			<< "Argument 1: " << "[path to the input graph]\n"
			<< "Argument 2 " << "[path to the  output image]\n";
	  return 1;
	}
	// Input parameters
	const std::string filename = argv[1];
	const std::string outputname = argv[2];

	using InputImageType               = otb::Image<float, 2>;
	using ImageToBaatz				   = otb::obia::ImageToBaatzGraphFilter<InputImageType>;
	using GraphType                    = ImageToBaatz::OutputGraphType;
	using GraphReader		  		   = otb::obia::GraphFileReader<GraphType>;
	using LabelPixelType               = unsigned int;
	using LabelImageType               = otb::Image< LabelPixelType, 2 >;
	using GraphToLabelImageFilterType  = otb::obia::GraphToLabelImageFilter<GraphType, LabelImageType>;
	using RGBPixelType                 = itk::RGBPixel<unsigned char>;
	using RGBImageType                 = otb::Image<RGBPixelType, 2>;
	using LabelToRGBFilterType         = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
	using RGBWriterType                = otb::ImageFileWriter< RGBImageType >;
	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;


	/* Pipeline setup */
	// Read the input image
	auto graphReader = GraphReader::New();
	graphReader->SetFileName(filename);
	// Convert contours to contour labels
	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	graphToLabelFilter->SetInput(graphReader->GetOutput());
	// Convert contour labels to image labels
	auto fillHoleFilter = FillholeFilterType::New();
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	// Convert to RGB image
	auto labelToRGBFilter = LabelToRGBFilterType::New();
	labelToRGBFilter->SetInput(fillHoleFilter->GetOutput());
	// Write RGB Image
	auto rgbWriter = RGBWriterType::New();
	rgbWriter->SetFileName(outputname);
	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
	rgbWriter->Update();
	return EXIT_SUCCESS;
}




