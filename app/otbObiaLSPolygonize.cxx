#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbObiaLSPolygonizeScheduler.h"
#include "otbObiaDouglasPeukerSimplify.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaConstExpr.h"
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

	std::cout<<"Init"<<std::endl;
        //General description
        SetName("LSPolygonize");
        SetDescription("Large Scale Image Polygonization Application");

        //Documentation
        SetDocName("Large Scale Polygonization");
        SetDocLongDescription("This application provides several methods to perform polygonization of very high resolution images");
        SetDocLimitations("None");
        SetDocAuthors("OBIA-Team");
        SetDocSeeAlso(" ");

	std::cout<<"IO"<<std::endl;
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


	// Processing parameters
	std::cout<<"Processing"<<std::endl;
        AddParameter(ParameterType_Group,"processing","Set of parameters related to processing");
	MandatoryOff("processing");
        AddParameter(ParameterType_Int,  "processing.nbtilesx",   "Number of tiles X");
        SetParameterDescription( "processing.nbtilesx", "Number of tiles X");
	SetDefaultParameterInt("processing.nbtilesx",1);
        AddParameter(ParameterType_Int,  "processing.nbtilesy",   "Number of tiles Y");
        SetParameterDescription( "processing.nbtilesy", "Number of tiles Y");
	SetDefaultParameterInt("processing.nbtilesy",1);

        AddParameter(ParameterType_Int,  "processing.maxtilesizex",   "Max tile size X");
        SetParameterDescription( "processing.maxtilesizex",   "Max tile size X");
	SetDefaultParameterInt("processing.maxtilesizex",250);
        AddParameter(ParameterType_Int,  "processing.maxtilesizey",   "Max tile size Y");
        SetParameterDescription( "processing.maxtilesizey",   "Max tile size Y");
	SetDefaultParameterInt("processing.maxtilesizey",250);

    }

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

        uint32_t nbX = GetParameterInt("processing.nbtilesx");
        uint32_t nbY = GetParameterInt("processing.nbtilesy");
        uint32_t maxSizeX = GetParameterInt("processing.maxtilesizex");
        uint32_t maxSizeY = GetParameterInt("processing.maxtilesizey");
	std::cout<<"Parametres recuperes"<<std::endl;
        // Polygonize
        using InputGraphType = otb::obia::Graph< otb::obia::Node<
                                                 otb::obia::BaatzNodeAttribute,
                                                 otb::obia::BaatzEdgeAttribute> >;
        using TilesMapType = std::map< int, std::set<uint32_t> >;
        using SimplifyFuncType =  otb::obia::DouglasPeukerFunc;
        using LSPolygonizeSchedulerType = otb::obia::LSPolygonizeScheduler<InputGraphType, SimplifyFuncType>;

	std::cout<<"Filtres declares"<<std::endl;
	std::cout<<"Graph to be used"<<filename<<std::endl;
        //Read graph from disk
        auto graph = otb::obia::GraphOperations<InputGraphType>::ReadGraphFromDisk(filename);

        //Create filter to vectorize
        auto LSPolygonizeScheduler = LSPolygonizeSchedulerType::New();

        auto simplifyFunc = new SimplifyFuncType;
        simplifyFunc->SetTolerance(1.0);

	size_t found;
	std::cout << "Splitting: " << filename << std::endl;
	found=filename.find_last_of("/\\");

        //Number of tiles per processor
        TilesMapType tilesPerProcessor;
        tilesPerProcessor[0].insert(0);
        tilesPerProcessor[1].insert(0);
	std::cout<<"Large scale processing parameters"<<std::endl;
        LSPolygonizeScheduler->SetImageHeight(graph->GetImageWidth());
        LSPolygonizeScheduler->SetImageWidth(graph->GetImageWidth());
        LSPolygonizeScheduler->SetNumberOfTilesX(nbX);
        LSPolygonizeScheduler->SetNumberOfTilesY(nbY);
        LSPolygonizeScheduler->SetMaxTileSizeX(maxSizeX);
        LSPolygonizeScheduler->SetMaxTileSizeY(maxSizeY);

        LSPolygonizeScheduler->SetTilesPerProcessor(tilesPerProcessor);
        LSPolygonizeScheduler->SetSimplifyFunc(simplifyFunc);
        LSPolygonizeScheduler->SetWriteVector(true);
	LSPolygonizeScheduler->SetInputDirectory(filename.substr(0,found));
        LSPolygonizeScheduler->SetOutputDir(outDir);
        LSPolygonizeScheduler->SetTemporaryDirectory(tmpDir);

        LSPolygonizeScheduler->SetGraphPrefixName("Graph_Vector");
        LSPolygonizeScheduler->SetGraph(graph);
	std::cout<<"Update"<<std::endl;
        LSPolygonizeScheduler->Update();

        std::cout << "End application" << std::endl;

        //Get output in order to write OGRDS into a file?

    }
};

OTB_APPLICATION_EXPORT(LSPolygonize)
}
}
