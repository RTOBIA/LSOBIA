#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbObiaGraphOperations.h"
#include "otbObiaGraphToVectorFilter.h"
#include "otbObiaSimplifyVectorFilter.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaDouglasPeukerSimplify.h"

#include <string>
#include <sstream>
#define NUM_ELEMENT 4

using std::string;
using std::stringstream;

namespace otb
{
namespace Wrapper
{

class LSPolygonize : public Application
{
public:

	typedef LSPolygonize Self;
        typedef Application SuperClass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(LSPolygonize, Application);

private:


	// Init App
	void DoInit()
	{

		//General description
		SetName("LSPolygonize");
		SetDescription("Large Scale Image Polygonize Application");

		//Documentation
		SetDocName("Large Scale Polygonization");
		SetDocLongDescription("This application provides several methods to perform polygonization of very high resolution images");
		SetDocLimitations("None");
		SetDocAuthors("OBIA-Team");
		SetDocSeeAlso(" ");

		// IO Parameters
		AddParameter(ParameterType_Group,"io","Set of parameters related to input/output");
		AddParameter(ParameterType_String,  "io.gr",   "Input graph path");
		SetParameterDescription("io.gr", "Graph");
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
		std::string filename = GetParameterString("io.gr");
		std::string outDir = GetParameterString("io.out.dir");
		std::string gmlFile = GetParameterString("io.out.gmlfile");
		std::string tmpDir = GetParameterString("io.temp");

		// Polygonize
		using InputGraphType = otb::obia::Graph< otb::obia::Node<
												 otb::obia::BaatzNodeAttribute,
												 otb::obia::BaatzEdgeAttribute> >;

		using GraphToVectorFilterType = otb::obia::GraphToVectorFilter<InputGraphType>;

		//Read graph from disk
		auto graph = otb::obia::GraphOperations<InputGraphType>::ReadGraphFromDisk(filename);

		//Create filter to vectorize
		auto graphToVectorFilter = GraphToVectorFilterType::New();
		graphToVectorFilter->SetInput(graph);


		//Create filter to simplify
		using SimplifyVectorFilterType = otb::obia::SimplifyVectorFilter<otb::obia::DouglasPeukerFunc>;
		auto simplifyVectorFilter = SimplifyVectorFilterType::New();
		auto simplifyFunc = new otb::obia::DouglasPeukerFunc();
		simplifyVectorFilter->SetSimplifyFunc(simplifyFunc);
		simplifyVectorFilter->SetInput(graphToVectorFilter->GetOutput());
		simplifyVectorFilter->Update();

		//graphToVectorFilter->Update();
		std::cout << "End application" << std::endl;

		//Get outut in order to write OGRDS into a file

	}
};

OTB_APPLICATION_EXPORT(LSPolygonize)
}
}
