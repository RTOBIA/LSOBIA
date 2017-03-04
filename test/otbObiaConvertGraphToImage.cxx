/*
 * otbObiaConvertGraphToImage.cxx
 *
 *  Created on: 9 févr. 2017
 *      Author: isnard
 */
#include "otbVectorImage.h"
#include "otbImage.h"
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaGraphOperations.h"
#include "otbObiaGraphToLabelImageFilter.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImageFileWriter.h"
#include "gdal_priv.h"
#include "gdal_utils.h"


/** Some convenient alias */
using InputImageType = otb::VectorImage<float, 2>;
using ImageToBaatz				   = otb::obia::ImageToBaatzGraphFilter<InputImageType>;
using GraphType                    = ImageToBaatz::OutputGraphType;
using GraphPointerType            = typename GraphType::Pointer;
using GraphOperationsType         = otb::obia::GraphOperations<GraphType>;


using LabelPixelType = unsigned int;
using LabelImageType = otb::Image< LabelPixelType, 2 >;
using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<GraphType, LabelImageType>;
using RGBPixelType = itk::RGBPixel<unsigned char>;
using RGBImageType = otb::Image<RGBPixelType, 2>;
using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
using RGBWriterType = otb::ImageFileWriter< RGBImageType >;

void RedefineColors(LabelToRGBFilterType* labelToRGB);

int main(int argc, char * argv[])
{
	if(argc < 2)
	{

	  std::cerr << "Usage " << argv[0] << ": \n"
			<< "Argument 1: " << "[path to the input data]\n"
			<< "Argument 2: " << "[path to the output image]\n" << std::endl;
	  return 1;
	}

	// Input parameters
	const std::string filename = argv[1];
	const std::string outputname  = argv[2];


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
	graphToLabelFilter->Update();

	std::cout<< "Graph to Label done" << std::endl;
	labelToRGBFilter->SetInput(graphToLabelFilter->GetOutput());
	RedefineColors(labelToRGBFilter);

	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
	rgbWriter->Update();

	// Read with GDAL
    GDALDataset  *poDataset;
    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpen( outputname.c_str(), GA_ReadOnly );
    if( poDataset == NULL )
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

    char* papszArgv[] = { "-of", "PNG", "-a_nodata", "0 0 0",  NULL };
    //char* papszArgv[] = { "-of", "PNG",  NULL };
    std::cout << "New options" <<std::endl;
    GDALTranslateOptions* psOptionsIn = GDALTranslateOptionsNew(papszArgv, NULL);
   // GDALTranslateOptions* psOptionsIn = GDALTranslateOptionsNew(papszArgv, NULL);

    std::cout << "Prepare to translate!" <<std::endl;
    GDALDatasetH labelDataset = GDALTranslate(pszDest.c_str(), poDataset, psOptionsIn, NULL);

    std::cout << "Translasted!" <<std::endl;
    //Free
   // GDALTranslateOptionsFree(psOptionsIn);
    GDALClose 	( labelDataset	) ;
    GDALClose 	( poDataset	) ;
}


void RedefineColors(LabelToRGBFilterType* labelToRGB)
{
	labelToRGB->ResetColors();
	labelToRGB->AddColor(255, 1, 1);
	labelToRGB->AddColor(1, 205, 1);
	labelToRGB->AddColor(1, 1, 255);
	labelToRGB->AddColor(1, 255, 255);
	labelToRGB->AddColor(255, 1, 255);
	labelToRGB->AddColor(255, 127, 1);
	labelToRGB->AddColor(1, 101, 1);
	labelToRGB->AddColor(138, 43, 226);
	labelToRGB->AddColor(139, 35, 35);
	labelToRGB->AddColor(1, 1, 128);
	labelToRGB->AddColor(139, 139, 1);
	labelToRGB->AddColor(255, 62, 151);
	labelToRGB->AddColor(139, 76, 57);
	labelToRGB->AddColor(1, 134, 139);
	labelToRGB->AddColor(205, 104, 57);
	labelToRGB->AddColor(191, 62, 255);
	labelToRGB->AddColor(1, 139, 69);
	labelToRGB->AddColor(199, 21, 133);
	labelToRGB->AddColor(205, 55, 1);
	labelToRGB->AddColor(32, 178, 171);
	labelToRGB->AddColor(106, 91, 205);
	labelToRGB->AddColor(255, 21, 147);
	labelToRGB->AddColor(69, 139, 116);
	labelToRGB->AddColor(72, 118, 255);
	labelToRGB->AddColor(205, 79, 57);
	labelToRGB->AddColor(1, 1, 205);
	labelToRGB->AddColor(139, 34, 82);
	labelToRGB->AddColor(139, 1, 139);
	labelToRGB->AddColor(238, 131, 238);
	labelToRGB->AddColor(139, 1, 1);
}


