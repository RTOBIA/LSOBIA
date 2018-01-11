#include "otbVectorImage.h"
#include "otbImage.h"
#include "otbObiaImageToBaatzGraphFilter.txx"
#include "otbObiaGraphOperations.txx"
#include "otbObiaGraphToLabelImageFilter.txx"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImageFileWriter.h"
#include "gdal_priv.h"
#include "gdal_utils.h"


int otbObiaConvertGraphToImage(int argc, char * argv[])
{

    /** Some convenient alias */
    using InputImageType = otb::VectorImage<float, 2>;
    using ImageToBaatz                   = otb::obia::ImageToBaatzGraphFilter<InputImageType>;
    using GraphType                    = ImageToBaatz::OutputGraphType;
    using GraphOperationsType         = otb::obia::GraphOperations<GraphType>;


    using LabelPixelType = unsigned int;
    using LabelImageType = otb::Image< LabelPixelType, 2 >;
    using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<GraphType, LabelImageType>;
    using RGBPixelType = itk::RGBPixel<unsigned char>;
    using RGBImageType = otb::Image<RGBPixelType, 2>;
    using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
    using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
    using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;


    if(argc < 2)
    {
        std::cout<<"Arg < 2..."<<std::endl;
        std::cerr << "Usage " << argv[0] << ": \n"
            << "Argument 1: " << "[path to the input data]\n"
            << "Argument 2: " << "[path to the output image]\n" << std::endl;
        return 1;
    }

    // Input parameters
    const std::string filename = argv[1];
    const std::string outputname  = argv[2];
    std::cout<<"Input parameters graph : "<<filename<<std::endl;
    std::cout<<"Output parameters image : "<<outputname<<std::endl;


    //Lecture graph
    std::cout << "Lecture Graphe : " << filename << std::endl;
    auto graph = GraphOperationsType::ReadGraphFromDisk(filename);
    graph->SetImageHeight(1000);
    graph->SetImageWidth(1000);
    graph->SetNumberOfSpectralBands(4);
    std::cout << "Nombre de noeud : " << graph->GetNumberOfNodes() << std::endl;
    std::cout << "Ecriture vers : " << outputname << std::endl;


    auto graphToLabelFilter = GraphToLabelImageFilterType::New();
    auto labelToRGBFilter = LabelToRGBFilterType::New();
    auto rgbWriter = RGBWriterType::New();
    rgbWriter->SetFileName(outputname);
    graphToLabelFilter->SetInput(graph);
    //graphToLabelFilter->Update();
    
    auto fillHoleFilter = FillholeFilterType::New();
    fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());

    std::cout<< "Graph to Label done" << std::endl;
    labelToRGBFilter->SetInput(fillHoleFilter->GetOutput());


    labelToRGBFilter->ResetColors();
    labelToRGBFilter->AddColor(255, 1, 1);
    labelToRGBFilter->AddColor(1, 205, 1);
    labelToRGBFilter->AddColor(1, 1, 255);
    labelToRGBFilter->AddColor(1, 255, 255);
    labelToRGBFilter->AddColor(255, 1, 255);
    labelToRGBFilter->AddColor(255, 127, 1);
    labelToRGBFilter->AddColor(1, 101, 1);
    labelToRGBFilter->AddColor(138, 43, 226);
    labelToRGBFilter->AddColor(139, 35, 35);
    labelToRGBFilter->AddColor(1, 1, 128);
    labelToRGBFilter->AddColor(139, 139, 1);
    labelToRGBFilter->AddColor(255, 62, 151);
    labelToRGBFilter->AddColor(139, 76, 57);
    labelToRGBFilter->AddColor(1, 134, 139);
    labelToRGBFilter->AddColor(205, 104, 57);
    labelToRGBFilter->AddColor(191, 62, 255);
    labelToRGBFilter->AddColor(1, 139, 69);
    labelToRGBFilter->AddColor(199, 21, 133);
    labelToRGBFilter->AddColor(205, 55, 1);
    labelToRGBFilter->AddColor(32, 178, 171);
    labelToRGBFilter->AddColor(106, 91, 205);
    labelToRGBFilter->AddColor(255, 21, 147);
    labelToRGBFilter->AddColor(69, 139, 116);
    labelToRGBFilter->AddColor(72, 118, 255);
    labelToRGBFilter->AddColor(205, 79, 57);
    labelToRGBFilter->AddColor(1, 1, 205);
    labelToRGBFilter->AddColor(139, 34, 82);
    labelToRGBFilter->AddColor(139, 1, 139);
    labelToRGBFilter->AddColor(238, 131, 238);
    labelToRGBFilter->AddColor(139, 1, 1);

    rgbWriter->SetInput(labelToRGBFilter->GetOutput());
    rgbWriter->Update();

    // Read with GDAL
    GDALDataset  *poDataset;
    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpen( outputname.c_str(), GA_ReadOnly );
    if( poDataset == nullptr )
    {
        std::cerr << "Erreur: impossible d'ouvrir " << outputname << std::endl;
    }

    //Convert all zeros value to no data
    std::string pszDest = outputname + "_no_data.png";

    //Setting options
   // char** papszArgv =   {"-a_nodata value 0"};
   /* char* papszArgv[] =
                            {
                                "-of", "tif",
                                "-a_nodata", "0",
                                // "-a_srs", "\"+proj=stere +lat_0=90 +lon_0=0 +lat_ts=60 +a=6378.14 +b=6356.75 +x_0=0 y_0=0\"",
                                //"-a_ullr", "0.0, -3650.000, 800.000, -4415.000"
                            };*/

    const char* papszArgv[] = { "-of", "PNG", "-a_nodata", "0 0 0",  nullptr };
    //char* papszArgv[] = { "-of", "PNG",  NULL };
    std::cout << "New options" <<std::endl;
    GDALTranslateOptions* psOptionsIn = GDALTranslateOptionsNew(const_cast<char**>(papszArgv), nullptr);
    // GDALTranslateOptions* psOptionsIn = GDALTranslateOptionsNew(papszArgv, NULL);

    std::cout << "Prepare to translate!" <<std::endl;
    GDALDatasetH labelDataset = GDALTranslate(pszDest.c_str(), poDataset, psOptionsIn, nullptr);

    std::cout << "Translasted!" <<std::endl;
    //Free
    // GDALTranslateOptionsFree(psOptionsIn);
    GDALClose     ( labelDataset    ) ;
    GDALClose     ( poDataset    ) ;

    return EXIT_SUCCESS;
}

