#include <otbObiaComputeAttributesFilter.h>
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbObiaGraphOperations.h"
#include "otbObiaGraphToVectorFilter.h"
#include "otbObiaConstExpr.h"
#include "otbImageFileReader.h"
#include <string>
#include <sstream>

#define NUM_ELEMENT 4

using std::string;
using std::stringstream;

namespace otb
{
namespace Wrapper
{

class GraphPolygonize : public Application
{
public:

	typedef GraphPolygonize Self;
        typedef Application SuperClass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(LSPolygonize, Application);

private:


	// Init App
	void DoInit()
	{

		//General description
		SetName("ComputeAttributs");
		SetDescription("Compute vector attributs Application");

		//Documentation
		SetDocName("Compute attributs");
		SetDocLongDescription("This application computes several attributs of a vector file");
		SetDocLimitations("None");
		SetDocAuthors("OBIA-Team");
		SetDocSeeAlso(" ");

		// IO Parameters
		AddParameter(ParameterType_Group,"io","Set of parameters related to input/output");
		AddParameter(ParameterType_String,  "io.vec",   "Input vector");
		SetParameterDescription("io.vec", "Vector");
		AddParameter(ParameterType_String,  "io.im",   "Input image");
		SetParameterDescription("io.im", "Image");
		AddParameter(ParameterType_Group, "io.out",  "Output directory");
		AddParameter(ParameterType_Directory, "io.out.dir",  "Output directory");
		SetParameterDescription("io.out.dir", "Output Directory");
		AddParameter(ParameterType_String, "io.out.gmlfile",  "GML FileName");
		SetParameterDescription("io.out.gmlfile", "GML FileName");
        MandatoryOff("io.out.gmlfile");
		AddParameter(ParameterType_Directory, "io.temp",  "Directory used for temporary data");
		SetParameterDescription("io.temp", "Temporary directory");



		/* TODO : remove this when the default values and choices have been implemented
		MandatoryOff("fusion.sylvester.linearcombination.image");
		AddParameter(ParameterType_Float,"fusion.glp.ratio","Resolutions ratio between the Panchromatic and the multispectral inputs");
		SetDefaultParameterFloat("",  4.);
		SetMinimumParameterFloatValue("", 0);
		SetDocExampleParameterValue("boolean", "true");
		SetDocExampleParameterValue("in", "QB_Suburb.png");
		SetDocExampleParameterValue("out", "Application_Example.png");
		*/
	}

	// TODO : parameter update should go there
	void DoUpdateParameters()
	{
	}

	// Execute App
	void DoExecute()
	{

		/* Global parameters */
		std::string vectorFile = GetParameterString("io.vec");
		std::string imagePath = GetParameterString("io.im");
		std::string outDir = GetParameterString("io.out.dir");
		std::string gmlFile = GetParameterString("io.out.gmlfile");
		std::string tmpDir = GetParameterString("io.temp");

		//Read image
		//Image reader
		using InputImageType              = otb::VectorImage<float, 2>;
		auto imgReader = otb::ImageFileReader<InputImageType>::New();
		imgReader->SetFileName(imagePath);

		//Compute attributs
		using ComputeAttributsFilterType = otb::obia::ComputeAttributesFilter<InputImageType>;
		auto computeAttributs = ComputeAttributsFilterType::New();

		const ogr::DataSource* outputDs =ogr::DataSource::New(vectorFile); //default 2nd argument is read
		computeAttributs->SetOGRData(outputDs);
		computeAttributs->SetInput(imgReader->GetOutput());
		computeAttributs->Update();
		//Read vector file

		//Create filter


//		std::cout << "Nombre layer = " << graphToVectorFilter->GetOutput()->GetLayersCount() << std::endl;
//		otb::ogr::Layer layer = graphToVectorFilter->GetOutput()->GetLayer(otb::obia::cleanedLayerName);
//		std::cout << "Get layer  cleaned = " << layer.GetFeatureCount(true) << std::endl;
		//graphToVectorFilter->Update();
		std::cout << "End application" << std::endl;

		//Get outut in order to write OGRDS into a file

	}
};

OTB_APPLICATION_EXPORT(GraphPolygonize)
}
}
