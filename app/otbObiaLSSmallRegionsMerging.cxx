#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

#include "otbObiaLSBaatzSegmentationScheduler.h"
#include "otbObiaLSSmallRegionsMergingScheduler.h"
#include "otbObiaLSMeanShiftScheduler.h"

#include "otbVectorImage.h"

#include <string>
#include <sstream>
#define NUM_ELEMENT 4

using std::string;
using std::stringstream;

namespace otb
{
namespace Wrapper
{

class LSSmallRegionsMerging : public Application
{
public:

	typedef LSSmallRegionsMerging Self;
        typedef Application SuperClass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(LSSmallRegionsMerging, Application);

private:

	// Available algorithms
	enum Algorithm
	{
		ALG_BAATZ,
		ALG_MEANSHIFT
	};

	// Available search modes
	enum SearchMode
	{
		SM_ON,
		SM_OFF
	};

	// Init App
	void DoInit()
	{

		//General description
		SetName("LSSmallRegionsMerging");
		SetDescription("Large Scale Image Small Regions Merging Application");

		//Documentation
		SetDocName("Large Scale Small Regions Merging");
		SetDocLongDescription("This application provides a method to perform small regions merging of very high resolution images");
		SetDocLimitations("None");
		SetDocAuthors("OBIA-Team");
		SetDocSeeAlso(" ");

		// IO Parameters
		AddParameter(ParameterType_Group,"io","Set of parameters related to input/output");
		AddParameter(ParameterType_String,  "io.im",   "Input image path");
		SetParameterDescription("io.im", "Image");
		AddParameter(ParameterType_Group, "io.out",  "Output directory");
		AddParameter(ParameterType_Directory, "io.out.dir",  "Output directory");
		SetParameterDescription("io.out.dir", "Output Directory");
		AddParameter(ParameterType_String, "io.out.labelimage",  "Label Image Name");
		SetParameterDescription("io.out.labelimage", "Label Image Name");
                MandatoryOff("io.out.labelimage");
		AddParameter(ParameterType_Directory, "io.temp",  "Directory used for temporary data");
		SetParameterDescription("io.temp", "Temporary directory");

		// Algorithm Parameters
		AddParameter(ParameterType_Choice,"algorithm","Segmentation algorithm name");
		AddChoice("algorithm.baatz", "Baatz and Shape algorithm" );
		AddChoice("algorithm.meanshift", "Mean-shift algorithm" );

		AddParameter(ParameterType_Int,"algorithm.baatz.numitfirstpartial","Number of iterations for first partial segmentation");
		AddParameter(ParameterType_Int,"algorithm.baatz.numitpartial","Number of iterations for partial segmentation");
		AddParameter(ParameterType_Float,"algorithm.baatz.stopping","Value for stopping criterion");
		AddParameter(ParameterType_Float,"algorithm.baatz.spectralweight","Value for spectral weight");
		AddParameter(ParameterType_Float,"algorithm.baatz.geomweight","Value for geometric (shape) weight");
		//AddParameter(ParameterType_InputVectorData,"baatz.bandweights", "optional band weights");
		//MandatoryOff("baatz.bandweights");

		AddParameter(ParameterType_Int,"algorithm.meanshift.maxiter","max number of iterations");
		AddParameter(ParameterType_Float,"algorithm.meanshift.spatialr","Spatial bandwidth");
		AddParameter(ParameterType_Float,"algorithm.meanshift.spectralr","Spectral bandwidth");
		AddParameter(ParameterType_Float,"algorithm.meanshift.threshold","Threshold");
		AddParameter(ParameterType_Float,"algorithm.meanshift.ranger","Spectral range ramp");
		AddParameter(ParameterType_Choice,"algorithm.meanshift.modesearch","Activation of search mode");
		AddChoice("algorithm.meanshift.modesearch.on","Activated");
		AddChoice("algorithm.meanshift.modesearch.off","Deactivated");

		// Processing Parameters
		AddParameter(ParameterType_Group,"processing","Set of parameters related to parallel processing options");
		AddParameter(ParameterType_Int,"processing.memory","Maximum memory to be used on the main node");
		AddParameter(ParameterType_Int,"processing.maxtilesizex","Maximum size of tiles along x axis");
		AddParameter(ParameterType_Int,"processing.maxtilesizey","Maximum size of tiles along x axis");
		AddParameter(ParameterType_Int,"processing.minimalsurface","Minimal surface");

	}

	// TODO : parameter update should go there
	void DoUpdateParameters()
	{
	}

	// Execute App
	void DoExecute()
	{
		/* Global parameters */

		std::string filename = GetParameterString("io.im");
		std::string outDir = GetParameterString("io.out.dir");
		std::string labelImage = GetParameterString("io.out.labelimage");
		std::string tmpDir = GetParameterString("io.temp");

		uint32_t maxTileWidth = GetParameterInt("processing.maxtilesizex");
		uint32_t maxTileHeight = GetParameterInt("processing.maxtilesizey");
		unsigned long int memory = GetParameterInt("processing.memory");



		// TODO : factorize parameter setup to avoid creating both filters
		// BAATZ
		using InputImageType = otb::VectorImage<float, 2>;
		using LSBaatzSegmentationSchedulerType = otb::obia::LSBaatzSegmentationScheduler<InputImageType>;
		auto lsBaatzFilter = LSBaatzSegmentationSchedulerType::New();
		uint32_t nbStartingIterations;
		uint32_t nbPartialIterations;
		float threshold;
		float spectralW;
		float shapeW;
		std::vector<float> bandWeights;

		// MEANSHIFT
		using LabelPixelType              = unsigned int;
		using LSMeanShiftSchedulerType       = otb::obia::LSMeanShiftScheduler<InputImageType, LabelPixelType>;
		auto lsMSFilter = LSMeanShiftSchedulerType::New();
		unsigned int maxIter;
		unsigned int spatialr;
		float spectralr;
		float ranger;
		bool modeSearch = false;

		//SRM
		using LSSmallRegionsMergingType = otb::obia::LSSmallRegionsMergingScheduler<LSBaatzSegmentationSchedulerType::OutputGraphType>;
		auto lsSMR = LSSmallRegionsMergingType::New();

		switch(GetParameterInt("algorithm"))
		{
			case ALG_BAATZ :
			{
				nbStartingIterations = GetParameterInt("algorithm.baatz.numitfirstpartial");
				nbPartialIterations = GetParameterInt("algorithm.baatz.numitpartial");
				threshold = GetParameterFloat("algorithm.baatz.stopping");
				threshold = threshold * threshold;
				spectralW = GetParameterFloat("algorithm.baatz.spectralweight");
				shapeW = GetParameterFloat("algorithm.baatz.geomweight");

				lsBaatzFilter->SetFileName(filename);
				lsBaatzFilter->SetMaxTileSizeX(maxTileWidth);
				lsBaatzFilter->SetMaxTileSizeY(maxTileHeight);
				lsBaatzFilter->SetStartingNumberOfIterations(nbStartingIterations);
				lsBaatzFilter->SetPartialNumberOfIterations(nbPartialIterations);
				lsBaatzFilter->SetTemporaryDirectory(tmpDir);
				lsBaatzFilter->SetAvailableMemory(memory);
				lsBaatzFilter->SetThreshold(threshold);
				lsBaatzFilter->SetSpectralWeight(spectralW);
				lsBaatzFilter->SetShapeWeight(shapeW);
				lsBaatzFilter->SetWriteLabelImage(true);
				lsBaatzFilter->SetWriteGraph(true);
                lsBaatzFilter->SetLabelImageName(labelImage);
				lsBaatzFilter->SetOutputDir(outDir);
				lsBaatzFilter->Update();

				//Post processing

				unsigned long int minimalSurface = GetParameterInt("processing.minimalsurface");
				lsSMR->SetAvailableMemory(memory);
				lsSMR->SetOutputDir(outDir);
				lsSMR->SetTemporaryDirectory(tmpDir);
				lsSMR->SetMaxTileSizeX(maxTileWidth);
				lsSMR->SetMaxTileSizeY(maxTileHeight);


				lsSMR->SetImageHeight(lsBaatzFilter->GetImageHeight());
				lsSMR->SetImageWidth(lsBaatzFilter->GetImageWidth());
				lsSMR->SetNumberOfSpectralBands(lsBaatzFilter->GetNumberOfSpectralBands());
				lsSMR->SetTileMap(lsBaatzFilter->GetTileMap());
				lsSMR->SetTilesPerProcessor(lsBaatzFilter->GetTilesPerProcessor());
				lsSMR->SetNumberOfTilesX(lsBaatzFilter->GetNumberOfTilesX());
				lsSMR->SetNumberOfTilesY(lsBaatzFilter->GetNumberOfTilesY());


				lsSMR->SetMaxNumberOfTilesPerProcessor(lsBaatzFilter->GetMaxNumberOfTilesPerProcessor());
				lsSMR->SetMinimalSurface(minimalSurface);
				lsSMR->SetNumberOfIterations(nbPartialIterations);
				lsSMR->SetAggregateGraph(false);
				lsSMR->Update();

				break;
			}
			case ALG_MEANSHIFT :
			{
				maxIter = GetParameterInt("algorithm.meanshift.maxiter");
				spatialr = GetParameterFloat("algorithm.meanshift.spatialr");
				spectralr = GetParameterFloat("algorithm.meanshift.spectralr");
				threshold = GetParameterFloat("algorithm.meanshift.threshold");
				ranger = GetParameterFloat("algorithm.meanshift.ranger");
				switch (GetParameterInt("algorithm.meanshift.modesearch")) {
					case SM_ON:
						modeSearch = true;
						break;
					case SM_OFF:
						modeSearch = false;
						break;
				}

				lsMSFilter->SetFileName(filename);
				lsMSFilter->SetMaxTileSizeX(maxTileWidth);
				lsMSFilter->SetMaxTileSizeY(maxTileHeight);
				lsMSFilter->SetAvailableMemory(memory);
				lsMSFilter->SetTemporaryDirectory(tmpDir);
				lsMSFilter->SetMaxNumberOfIterations(maxIter);
				lsMSFilter->SetSpatialBandWidth(spatialr);
				lsMSFilter->SetSpectralRangeBandWidth(spectralr);
				lsMSFilter->SetThreshold(threshold);
				lsMSFilter->SetSpectralRangeRamp(ranger);
				lsMSFilter->SetModeSearch(modeSearch);
				lsMSFilter->SetOutputDir(outDir);
                lsMSFilter->SetLabelImageName(labelImage);
				lsMSFilter->SetWriteLabelImage(true);
				lsMSFilter->SetWriteGraph(true);
				lsMSFilter->Update();
				break;
			}
			default:
			{
				break;
			}
		}

	}
};

OTB_APPLICATION_EXPORT(LSSmallRegionsMerging)
}
}
