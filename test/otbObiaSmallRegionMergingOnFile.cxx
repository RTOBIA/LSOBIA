#include <iostream>

#include "otbObiaBaatzSegmentationFilter.h"
#include "otbObiaBaatzGraphToSRMGraphFilter.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaLSMeanShiftScheduler.h"
#include "otbObiaSmallRegionsMergingFilter.h"



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
	using InputGraphType			  = otb::obia::Graph<otb::obia::Node < otb::obia::BaatzNodeAttribute,
																		 	 otb::obia::BaatzEdgeAttribute> >;
	using GraphPointerType 			  = typename InputGraphType::Pointer;
	using LabelPixelType              = unsigned int;
	using SmallRegionMergingFilter	= otb::obia::SmallRegionsMergingFilter<InputGraphType>;

	//Lecture gu graphe
	using GraphOperationsType         = otb::obia::GraphOperations<InputGraphType>;

	std::stringstream os;
	os << tmpDir << filename;
	auto graph = GraphOperationsType::ReadGraphFromDisk(os.str());

	std::cout << "Nombre noeud : " << graph->GetNumberOfNodes() << std::endl;

	std::cout << "Noeud 1 " << graph->GetNodeAt(0)->m_Attributes.m_Area << std::endl;

	auto SRMFilter  = SmallRegionMergingFilter::New();


	SRMFilter->SetInput(GraphPointerType (graph) );
	SRMFilter->SetMinimalSurface(minimalSurface);
	SRMFilter->UpdateLargestPossibleRegion();


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
	labelToRGBFilter->SetInput(graphToLabelFilter->GetOutput());
	//labelToRGBFilter->UpdateLargestPossibleRegion();
	rgbWriter->SetFileName(tmpDir + "/smallRegion.tif");
	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
	rgbWriter->Update();

	std::cout << "End Write" <<std::endl;

	  return EXIT_SUCCESS;
}

