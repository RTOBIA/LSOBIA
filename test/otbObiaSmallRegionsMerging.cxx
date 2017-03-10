#include <iostream>

#include "otbObiaBaatzSegmentationFilter.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaLSMeanShiftScheduler.h"
#include "otbObiaSmallRegionsMergingFilter.h"
#include "otbVectorImage.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbObiaGenericRegionMergingFilter.h"

int otbObiaSmallRegionsMerging(int argc, char * argv[])
{
	std::cout << "--------- Baatz + Small Region Merging -------------" << std::endl;

	if(argc < 12)
	{

		std::cerr << "Usage " << argv[0] << ": \n"
				<< "Argument 1: " << "[path to the input image]\n"
				<< "Argument 2: " << "[tile size with respect to x axis]\n"
				<< "Argument 3: " << "[tile size with respect to y axis]\n"
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
	// Input parameters
	const std::string filename = argv[1];
	const uint32_t tileWidth = atoi(argv[2]);
	const uint32_t tileHeight = atoi(argv[3]);
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

	using InputImageType              = otb::VectorImage<float, 2>;
	using InputGraphType				= otb::obia::Graph<otb::obia::Node < otb::obia::BaatzNodeAttribute,
			otb::obia::BaatzEdgeAttribute> >;
	using LabelPixelType              = unsigned int;

	using ImageToBaatz				= otb::obia::ImageToBaatzGraphFilter<InputImageType>;
	using BaatzSegmentation			= otb::obia::GenericRegionMergingFilter<InputGraphType, InputGraphType,
																			otb::obia::BaatzMergingCost<float, InputGraphType> ,
																			otb::obia::BaatzHeuristic<InputGraphType>,
																			otb::obia::BaatzUpdateAttribute<InputGraphType> >;
	using SmallRegionMergingFilter	= otb::obia::GenericRegionMergingFilter<InputGraphType, InputGraphType,
			
	using LabelImageType               = otb::Image< LabelPixelType, 2 >;
	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;

	auto imageToBaaz = ImageToBaatz::New();
	auto baatzFilter = BaatzSegmentation::New();
	auto SRMFilter  = SmallRegionMergingFilter::New();


	//Load image as a graph
	// Read the input image
	auto imgReader = otb::ImageFileReader<InputImageType>::New();
	imgReader->SetFileName(filename);

	// Segmentation filter
	baatzFilter->SetMaxNumberOfIterations(75);
	baatzFilter->GetMergingCostFunc()->SetThreshold(threshold);
	baatzFilter->GetMergingCostFunc()->SetSpectralWeight(spectralW);
	baatzFilter->GetMergingCostFunc()->SetShapeWeight(shapeW);
	baatzFilter->GetMergingCostFunc()->SetBandWeights(bandWeights);
	baatzFilter->GetHeuristicFunc()->SetThreshold(threshold);


	baatzFilter->SetThreshold(threshold);
	baatzFilter->SetSpectralWeight(spectralW);
	baatzFilter->SetShapeWeight(shapeW);
	baatzFilter->SetBandWeights(bandWeights);

	// Pipeline branching
	imageToBaaz->SetInput(imgReader->GetOutput());
	baatzFilter->SetInput(const_cast< InputGraphType * > (imageToBaaz->GetOutput() ));
	//baatzFilter->Update();

	//baatzToSRM->SetInput( baatzFilter->GetOutput());

	SRMFilter->SetInput(const_cast< InputGraphType * > (baatzFilter->GetOutput() ));
	SRMFilter->GetMergingCostFunc()->SetMinimalSurface(200);
	SRMFilter->SetMaxNumberOfIterations(2);

	SRMFilter->UpdateLargestPossibleRegion();
	//SRMFilter->Update();


	std::cout << "Nombre noeud : " << SRMFilter->GetOutput()->GetNumberOfNodes() << std::endl;
	using LabelPixelType1 = unsigned int;
	using LabelImageType = otb::Image< LabelPixelType1, 2 >;
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<InputGraphType, LabelImageType>;

	using RGBPixelType = itk::RGBPixel<unsigned char>;
	using RGBImageType = otb::Image<RGBPixelType, 2>;
	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;

	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	auto labelToRGBFilter = LabelToRGBFilterType::New();
	auto rgbWriter = RGBWriterType::New();

	graphToLabelFilter->SetInput(const_cast< InputGraphType * > (SRMFilter->GetOutput()));
	graphToLabelFilter->Update();
	auto fillHoleFilter = FillholeFilterType::New();
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	labelToRGBFilter->SetInput(fillHoleFilter->GetOutput());

	rgbWriter->SetFileName(tmpDir + "/smallRegion.tif");
	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
	rgbWriter->Update();

	std::cout << "End Write" <<std::endl;

	return EXIT_SUCCESS;
}
