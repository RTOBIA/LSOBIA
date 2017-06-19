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
std::vector<otb::obia::GenericAttribute<InputImageType>*> CreateAttributs();

bool isOn(const std::string & str)
{
	if (str == std::string("on"))
	{
		return true;
	}
	return false;
}

int otbLSOBIAComputeAttributsTest(int argc, char *argv[])
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
	imgReader->Update();//TODO/ Fix the fact to do update twice...

	//Compute attributs
	using OGRDataSourceType = otb::ogr::DataSource;
	using ComputeAttributesFilterType = otb::obia::ComputeAttributesFilter<InputImageType>;
	auto computeAttributs = ComputeAttributesFilterType::New();

	OGRDataSourceType::Pointer outputDs = otb::ogr::DataSource::New(vectorFile, OGRDataSourceType::Modes::Read); //default 2nd argument is read
	std::cout << "Number of layers = " << outputDs->GetLayersCount() << std::endl;
	std::cout << "Layer = " << outputDs->GetLayer(0).GetName() << std::endl;
	computeAttributs->SetOGRData(outputDs);
	computeAttributs->SetInput(imgReader->GetOutput());
	computeAttributs->SetInputLayerName(otb::obia::reconstructedLayerName);
	computeAttributs->SetLayerIndex(0);
	computeAttributs->SetOutLayerName(otb::obia::attributsLayerName);
	computeAttributs->SetFieldName(otb::obia::startingCoordsFieldName);
	computeAttributs->SetAttributs(CreateAttributs());
	computeAttributs->SetOutputDir(outDir);
	computeAttributs->SetOutputFilename(gmlFile);
	computeAttributs->Update();

    std::cout << "SUCCESS" << std::endl;
    return EXIT_SUCCESS;
}

std::vector<otb::obia::GenericAttribute<InputImageType>*> CreateAttributs()
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
