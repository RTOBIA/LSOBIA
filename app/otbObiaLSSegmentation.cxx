#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

#include "otbObiaLSBaatzSegmentationScheduler.h"
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

class LSSegmentation : public Application
{
public:

	typedef LSSegmentation Self;
        typedef Application SuperClass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(LSSegmentation, Application);

private:
	
	// Available algorithms
	enum Algorithm
	{
		ALG_BAATZ,
		ALG_MEANSHIFT
	};

	// Available search modes
	enum Modes
	{
		ON,
		OFF
	};
	
	// Init App
	void DoInit()
	{
	
		//General description
		SetName("LSSegmentation");
		SetDescription("Large Scale Image Segmentation Application");

		//Documentation
		SetDocName("Large Scale Segmentation");
		SetDocLongDescription("This application provides several methods to perform segmentation of very high resolution images");
		SetDocLimitations("None");
		SetDocAuthors("OBIA-Team");
		SetDocSeeAlso(" ");
		
		// IO Parameters
		AddParameter(ParameterType_Group,"io","Set of parameters related to input/output");
		AddParameter(ParameterType_String,  "io.im",   "Input image path");
		SetParameterDescription("io.im", "Image");
		AddParameter(ParameterType_Directory, "io.out",  "Output directory");
		SetParameterDescription("io.out", "Output Directory");
		AddParameter(ParameterType_Directory, "io.temp",  "Directory used for temporary data");
		SetParameterDescription("io.temp", "Temporary directory");

		// Algorithm Parameters
		AddParameter(ParameterType_Choice,"algorithm","Segmentation algorithm name");
		AddChoice("algorithm.baatz", "Baatz and Shape algorithm" );
		AddChoice("algorithm.meanshift", "Mean-shift algorithm" );

		AddParameter(ParameterType_Int,"algorithm.baatz.numitfirstpartial","Number of iterations for first partial segmentation");
		SetDefaultParameterInt("algorithm.baatz.numitfirstpartial",  1);
		AddParameter(ParameterType_Int,"algorithm.baatz.numitpartial","Number of iterations for partial segmentation");
		SetDefaultParameterInt("algorithm.baatz.numitpartial",  1);
		AddParameter(ParameterType_Float,"algorithm.baatz.stopping","Value for stopping criterion");
		SetDefaultParameterFloat("algorithm.baatz.stopping",  40.);
		AddParameter(ParameterType_Float,"algorithm.baatz.spectralweight","Value for spectral weight");
		SetDefaultParameterFloat("algorithm.baatz.spectralweight",  0.05);
		AddParameter(ParameterType_Float,"algorithm.baatz.geomweight","Value for geometric (shape) weight");
		SetDefaultParameterFloat("algorithm.baatz.geomweight",  0.95);

		//AddParameter(ParameterType_InputVectorData,"algorithm.baatz.bandweights", "optional band weights");
		//MandatoryOff("algorithm.baatz.bandweights");

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
		AddParameter(ParameterType_Choice,"processing.writeimages","Activation of image traces");
		AddChoice("processing.writeimages.on","Activated");
		AddChoice("processing.writeimages.off","Deactivated");
		AddParameter(ParameterType_Choice,"processing.writegraphs","Activation of graph traces");
		AddChoice("processing.writegraphs.on","Activated");
		AddChoice("processing.writegraphs.off","Deactivated");
		AddParameter(ParameterType_Choice,"processing.aggregategraphs","Aggregation of graph traces");
		AddChoice("processing.aggregategraphs.on","Activated");
		AddChoice("processing.aggregategraphs.off","Deactivated");


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
		std::string filename = GetParameterString("io.im");
		std::string outDir = GetParameterString("io.out");
		std::string tmpDir = GetParameterString("io.temp");
		
		uint32_t maxTileWidth = GetParameterInt("processing.maxtilesizex");
		uint32_t maxTileHeight = GetParameterInt("processing.maxtilesizey");
		unsigned long int memory = GetParameterInt("processing.memory");
		bool writeImages = false;
		bool writeGraphs = false;
		bool aggregateGraphs = false;

		switch (GetParameterInt("processing.writeimages")) {
			case ON:
				writeImages = true;
				break;
			case OFF:
				writeImages = false;
				break;
		}
		switch (GetParameterInt("processing.writegraphs")) {
			case ON:
				writeGraphs = true;
				break;
			case OFF:
				writeGraphs = false;
				break;
		}
		switch (GetParameterInt("processing.aggregategraphs")) {
			case ON:
				aggregateGraphs = true;
				break;
			case OFF:
				aggregateGraphs = false;
				break;
		}
		
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
				lsBaatzFilter->SetWriteLabelImage(writeImages);
				lsBaatzFilter->SetWriteGraph(writeGraphs);
				lsBaatzFilter->SetOutputDir(outDir);
				lsBaatzFilter->Update();

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
					case ON:
						modeSearch = true;
						break;
					case OFF:
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
				lsMSFilter->SetWriteLabelImage(writeImages);
				lsMSFilter->SetWriteGraph(writeGraphs);
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

OTB_APPLICATION_EXPORT(LSSegmentation)
}
}
