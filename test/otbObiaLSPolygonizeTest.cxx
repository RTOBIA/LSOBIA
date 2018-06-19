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
#include "otbObiaLSPolygonizeScheduler.txx"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaDouglasPeukerSimplify.h"
#include "otbVectorImage.h"

#define NUM_ELEMENT 4

bool isOn(const std::string & str)
{
	if (str == std::string("on"))
	{
		return true;
	}
	return false;
}

int otbObiaLSPolygonize(int argc, char *argv[])
{
    if(argc < 9)
    {

        std::cerr << "Usage " << argv[0] << ": \n"
            << "Argument 1: " << "[input graph directory]\n"
            << "Argument 2: " << "[temporary directory to store intermediate files for a node]\n"
            << "Argument 3: " << "[output directory path (must me already created)]\n"
			<< "Argument 4: " << "[Image width]"
			<< "Argument 5: " << "[Image height]"
			<< "Argument 6: " << "[Number of tiles X]"
			<< "Argument 7: " << "[Number of tiles Y]"
			<< "Argument 8: " << "[Max tile size X]"
			<< "Argument 9: " << "[Max tile size Y]"
			<< std::endl;

        return 1;
    }

    /* Initialize MPI */
    auto mpiConfig = otb::MPIConfig::Instance();
    mpiConfig->Init(argc, argv);

    // Input parameters
    const std::string inDir  = argv[1];
    const std::string tmpDir = argv[2];
    const std::string outDir = argv[3];
    uint32_t imageWidth 	 = atoi(argv[4]);
    uint32_t imageHeight	 = atoi(argv[5]);
    uint32_t nbX			 = atoi(argv[6]);
    uint32_t nbY 			 = atoi(argv[7]);
    uint32_t maxSizeX		 = atoi(argv[8]);
    uint32_t maxSizeY 		 = atoi(argv[9]);

    std::cout << "Input directory = " << inDir << std::endl;
	// Polygonize
	using InputGraphType = otb::obia::Graph< otb::obia::Node<
											 otb::obia::BaatzNodeAttribute,
											 otb::obia::BaatzEdgeAttribute> >;
	using TilesMapType = std::map< int, std::set<uint32_t> >;
	using SimplifyFuncType =  otb::obia::DouglasPeukerFunc;
	using LSPolygonizeSchedulerType = otb::obia::LSPolygonizeScheduler<InputGraphType, SimplifyFuncType>;

	//Read graph from disk
	//auto graph = otb::obia::GraphOperations<InputGraphType>::ReadGraphFromDisk(filename);

	//Create filter to vectorize
	auto LSPolygonizeScheduler = LSPolygonizeSchedulerType::New();

	auto simplifyFunc = new SimplifyFuncType;
	simplifyFunc->SetTolerance(2.0);

	//Number of tiles per processor
//	TilesMapType tilesPerProcessor;
//	tilesPerProcessor[0].insert(0);
//	tilesPerProcessor[1].insert(0);

    //LSPolygonizeScheduler->SetTilesPerProcessor(tilesPerProcessor);
	LSPolygonizeScheduler->SetImageHeight(imageWidth);
	LSPolygonizeScheduler->SetImageWidth(imageHeight);
	LSPolygonizeScheduler->SetNumberOfTilesX(nbX);
	LSPolygonizeScheduler->SetNumberOfTilesY(nbY);
	LSPolygonizeScheduler->SetMaxTileSizeX(maxSizeX);
	LSPolygonizeScheduler->SetMaxTileSizeY(maxSizeY);

    LSPolygonizeScheduler->SetSimplifyFunc(simplifyFunc);
    LSPolygonizeScheduler->SetWriteVector(true);
    LSPolygonizeScheduler->SetIsSimplify(true);
    LSPolygonizeScheduler->SetInputDirectory(inDir);
    LSPolygonizeScheduler->SetOutputDir(outDir);
    LSPolygonizeScheduler->SetTemporaryDirectory(tmpDir);
    LSPolygonizeScheduler->SetGraphPrefixName("Graph");
    LSPolygonizeScheduler->Update();


    std::cout << "SUCCESS" << std::endl;
    return EXIT_SUCCESS;
}
