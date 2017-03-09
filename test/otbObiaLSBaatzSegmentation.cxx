#include "otbObiaLSBaatzSegmentationScheduler.txx"
#include "otbVectorImage.h"

#define NUM_ELEMENT 4

int otbObiaLSBaatzSegmentation(int argc, char *argv[])
{
	if(argc < 12)
	{

	  std::cerr << "Usage " << argv[0] << ": \n"
		    << "Argument 1: " << "[path to the input image]\n"
		    << "Argument 2: " << "[max tile size with respect to x axis]\n"
		    << "Argument 3: " << "[max tile size with respect to y axis]\n"
		    << "Argument 4: " << "[number of iterations for the first partial segmentation]\n"
		    << "Argument 5: " << "[number of iterations for the partial segmentations]\n"
		    << "Argument 6: " << "[temporary directory to store intermediate files for a node]\n"
		    << "Argument 7: " << "[available memory in the master node in megabytes]\n"
		    << "Argument 8: " << "[stopping criterion threshold]\n"
		    << "Argument 9: " << "[spectral weight]\n"
		    << "Argument 10: " << "[shape weight]\n"
		    << "Argument 11: " << "[output directory path (must me already created)]\n"
		    << "Argument 12: " << "[weights for each spectral band (optional)...]\n";
	  return 1;
	}

	/* Initialize MPI */
	auto mpiConfig = otb::MPIConfig::Instance();
	mpiConfig->Init(argc, argv);

	using InputImageType = otb::VectorImage<float, 2>;
	using LSBaatzSegmentationSchedulerType = otb::obia::LSBaatzSegmentationScheduler<InputImageType>;
	auto lsBaatzFilter = LSBaatzSegmentationSchedulerType::New();

	// Input parameters
	const std::string filename = argv[1];
	const uint32_t MaxTileWidth = atoi(argv[2]);
	const uint32_t MaxTileHeight = atoi(argv[3]);
	const uint32_t nbStartingIterations = atoi(argv[4]);
	const uint32_t nbPartialIterations = atoi(argv[5]);
	const std::string tmpDir = argv[6];
	const unsigned long int mem = atoi(argv[7]);
	const float threshold = atof(argv[8]) * atof(argv[8]);
	const float spectralW = atof(argv[9]);
	const float shapeW = atof(argv[10]);
	const std::string outDir = argv[11];
	std::vector<float> bandWeights;

	if(argc >= 13)
	  {
	    for(int i = 12; i < argc; i++)
	      {
		bandWeights.push_back(atof(argv[i]));
	      }
	  }

	lsBaatzFilter->SetFileName(filename);
	lsBaatzFilter->SetMaxTileSizeX(MaxTileWidth);
	lsBaatzFilter->SetMaxTileSizeY(MaxTileHeight);
	lsBaatzFilter->SetStartingNumberOfIterations(nbStartingIterations);
	lsBaatzFilter->SetPartialNumberOfIterations(nbPartialIterations);
	lsBaatzFilter->SetTemporaryDirectory(tmpDir);
	lsBaatzFilter->SetAvailableMemory(mem);
	lsBaatzFilter->SetThreshold(threshold);
	lsBaatzFilter->SetSpectralWeight(spectralW);
	lsBaatzFilter->SetShapeWeight(shapeW);
	lsBaatzFilter->SetBandWeights(bandWeights);
	lsBaatzFilter->SetWriteLabelImage(true);
	lsBaatzFilter->SetWriteGraph(true);
	lsBaatzFilter->SetOutputDir(outDir);
	lsBaatzFilter->Update();

	std::cout << "SUCCESS" << std::endl;
	return EXIT_SUCCESS;
}
