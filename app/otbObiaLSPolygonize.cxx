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

        AddParameter(ParameterType_Int,  "io.im.width",   "Image width");
        SetParameterDescription( "io.im.width", "Image width");
        AddParameter(ParameterType_Int,  "io.im.height",   "Image height");
        SetParameterDescription( "io.im.height", "height width");

        AddParameter(ParameterType_Int,  "io.im.nbTilesX",   "Number of tiles X");
        SetParameterDescription( "io.im.nbTilesX", "Number of tiles X");
        AddParameter(ParameterType_Int,  "io.im.nbTilesY",   "Number of tiles Y");
        SetParameterDescription( "io.im.nbTilesY", "Number of tiles Y");

        AddParameter(ParameterType_Int,  "io.im.maxTileSizeX",   "Max tile size X");
        SetParameterDescription( "io.im.maxTileSizeX",   "Max tile size X");
        AddParameter(ParameterType_Int,  "io.im.maxTileSizeY",   "Max tile size Y");
        SetParameterDescription( "io.im.maxTileSizeY",   "Max tile size Y");

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

        uint32_t imageWidth = GetParameterInt("io.im.width");
        uint32_t imageHeight = GetParameterInt("io.im.height");
        uint32_t nbX = GetParameterInt("io.im.nbTilesX");
        uint32_t nbY = GetParameterInt("io.im.nbTilesY");
        uint32_t maxSizeX = GetParameterInt("io.im.maxTileSizeX");
        uint32_t maxSizeY = GetParameterInt("io.im.maxTileSizeY");

        // Polygonize
        using InputGraphType = otb::obia::Graph< otb::obia::Node<
                                                 otb::obia::BaatzNodeAttribute,
                                                 otb::obia::BaatzEdgeAttribute> >;
        using TilesMapType = std::map< int, std::set<uint32_t> >;
        using SimplifyFuncType =  otb::obia::DouglasPeukerFunc;
        using LSPolygonizeSchedulerType = otb::obia::LSPolygonizeScheduler<InputGraphType, SimplifyFuncType>;

        //Read graph from disk
        auto graph = otb::obia::GraphOperations<InputGraphType>::ReadGraphFromDisk(filename);

        //Create filter to vectorize
        auto LSPolygonizeScheduler = LSPolygonizeSchedulerType::New();

        auto simplifyFunc = new SimplifyFuncType;
        simplifyFunc->SetTolerance(1.0);

        //Number of tiles per processor
        TilesMapType tilesPerProcessor;
        tilesPerProcessor[0].insert(0);
        tilesPerProcessor[1].insert(0);

        LSPolygonizeScheduler->SetImageHeight(imageWidth);
        LSPolygonizeScheduler->SetImageWidth(imageHeight);
        LSPolygonizeScheduler->SetNumberOfTilesX(nbX);
        LSPolygonizeScheduler->SetNumberOfTilesY(nbY);
        LSPolygonizeScheduler->SetMaxTileSizeX(maxSizeX);
        LSPolygonizeScheduler->SetMaxTileSizeY(maxSizeY);

        LSPolygonizeScheduler->SetTilesPerProcessor(tilesPerProcessor);
        LSPolygonizeScheduler->SetSimplifyFunc(simplifyFunc);
        LSPolygonizeScheduler->SetWriteVector(true);
        LSPolygonizeScheduler->SetOutputDir(outDir);
        LSPolygonizeScheduler->SetTemporaryDirectory(tmpDir);

        LSPolygonizeScheduler->SetGraphPrefixName("Graph_Vector");
        LSPolygonizeScheduler->SetGraph(graph);
        LSPolygonizeScheduler->Update();

//        std::cout << "Nombre layer = " << graphToVectorFilter->GetOutput()->GetLayersCount() << std::endl;
//        otb::ogr::Layer layer = graphToVectorFilter->GetOutput()->GetLayer(otb::obia::cleanedLayerName);
//        std::cout << "Get layer  cleaned = " << layer.GetFeatureCount(true) << std::endl;
        //graphToVectorFilter->Update();
        std::cout << "End application" << std::endl;

        //Get outut in order to write OGRDS into a file

    }
};

OTB_APPLICATION_EXPORT(LSPolygonize)
}
}
