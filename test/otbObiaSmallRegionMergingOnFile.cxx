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
#include <iostream>

#include "otbObiaBaatzSegmentationFilter.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaLSMeanShiftScheduler.h"
#include "otbObiaSmallRegionsMergingFilter.h"
#include "otbObiaGenericRegionMergingFilter.h"

int otbObiaSmallRegionMergingOnFile(int argc, char * argv[])
{
    if(argc < 3)
    {

      std::cerr << "Usage " << argv[0] << ": \n"
            << "Argument 1: " << "[input graph]\n"
            << "Argument 2: " << "[Minimal Surface]\n"
            << "Argument 3: " << "[temporary directory to store intermediate files for a node]\n" << std::endl;
      return 1;
    }
    // Input parameters
    const std::string filename = argv[1];
    uint32_t minimalSurface = atoi(argv[2]);
    const std::string tmpDir = argv[3];


    using InputImageType              = otb::VectorImage<float, 2>;
    using InputGraphType              = otb::obia::Graph<otb::obia::Node < otb::obia::BaatzNodeAttribute,
                                                                              otb::obia::BaatzEdgeAttribute> >;
    using GraphPointerType               = typename InputGraphType::Pointer;
    using SmallRegionMergingFilter    = otb::obia::GenericRegionMergingFilter<InputGraphType, InputGraphType,
            otb::obia::SRMMergingCost<float, InputGraphType> ,
            otb::obia::SRMHeuristic<InputGraphType>,
            otb::obia::SRMUpdateAttribute<InputGraphType> >;
    using LabelPixelType = unsigned int;
    using LabelImageType = otb::Image< LabelPixelType, 2 >;
    using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<InputGraphType, LabelImageType>;

    using RGBPixelType = itk::RGBPixel<unsigned char>;
    using RGBImageType = otb::Image<RGBPixelType, 2>;
    using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
    using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
    using GraphOperationsType         = otb::obia::GraphOperations<InputGraphType>;


    auto graph = GraphOperationsType::ReadGraphFromDisk(filename);
    std::cout << "Nombre noeud : " << graph->GetNumberOfNodes() << std::endl;
    auto SRMFilter  = SmallRegionMergingFilter::New();
    std::cout << "SRM filter created"<< std::endl;
    SRMFilter->SetInput((GraphPointerType)graph);
    std::cout << "Graph input set"<< std::endl;
    SRMFilter->GetMergingCostFunc()->SetMinimalSurface(minimalSurface);
    std::cout << "Minimal Surface set : "<<minimalSurface<< std::endl;
    SRMFilter->UpdateLargestPossibleRegion();
    std::cout << "UpdateLargestPossibleRegion done"<< std::endl;

    auto graphToLabelFilter = GraphToLabelImageFilterType::New();
    auto labelToRGBFilter = LabelToRGBFilterType::New();
    auto rgbWriter = RGBWriterType::New();

    graphToLabelFilter->SetInput(const_cast< InputGraphType * > (SRMFilter->GetOutput()));
    graphToLabelFilter->Update();

    labelToRGBFilter->SetInput(graphToLabelFilter->GetOutput());
    rgbWriter->SetFileName(tmpDir + "/smallRegion.tif");
    rgbWriter->SetInput(labelToRGBFilter->GetOutput());
    rgbWriter->Update();

    std::cout << "End Write" <<std::endl;

    return EXIT_SUCCESS;
}

