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
    if(argc < 3)
    {

        std::cerr << "Usage " << argv[0] << ": \n"
            << "Argument 1: " << "[input graph directory]\n"
            << "Argument 2: " << "[temporary directory to store intermediate files for a node]\n"
            << "Argument 3: " << "[output directory path (must me already created)]\n" << std::endl;
        return 1;
    }

    /* Initialize MPI */
    auto mpiConfig = otb::MPIConfig::Instance();
    mpiConfig->Init(argc, argv);

    // Input parameters
    const std::string inDir  = argv[1];
    const std::string tmpDir = argv[2];
    const std::string outDir = argv[3];

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
