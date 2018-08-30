/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

#include "otbObiaLSBaatzSegmentationScheduler.h"
#include "otbObiaLSMeanShiftScheduler.h"

#include "otbImageFileReader.h"
#include "otbVectorImage.h"

#include "otbObiaSpatialTools.h"

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

	// Tiling configuration
    obia::TilingConfiguration tilingConf = {0, 0, 0};
    
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
        SetDefaultParameterInt("algorithm.baatz.numitfirstpartial",  1);
        MandatoryOff("algorithm.baatz.numitfirstpartial");
	
        AddParameter(ParameterType_Int,"algorithm.baatz.numitpartial","Number of iterations for partial segmentation");
        MandatoryOff("algorithm.baatz.numitpartial");
        SetDefaultParameterInt("algorithm.baatz.numitpartial",  1);

        AddParameter(ParameterType_Int,"algorithm.baatz.maxiter","max number of iterations");
        MandatoryOff("algorithm.baatz.maxiter");
        SetDefaultParameterInt("algorithm.baatz.maxiter",  75);

        AddParameter(ParameterType_Float,"algorithm.baatz.mindec","minimum decreasing of accumulated Memory");
        MandatoryOff("algorithm.baatz.mindec");
        SetDefaultParameterFloat("algorithm.baatz.mindec",  0.);

        AddParameter(ParameterType_Float,"algorithm.baatz.scale","Value for scale criterion");
        MandatoryOff("algorithm.baatz.scale");
        SetDefaultParameterFloat("algorithm.baatz.scale",  60.);
	
        AddParameter(ParameterType_Float,"algorithm.baatz.spectralweight","Value for spectral weight");
        SetDefaultParameterFloat("algorithm.baatz.spectralweight",  0.05);
        MandatoryOff("algorithm.baatz.spectralweight");

        AddParameter(ParameterType_Float,"algorithm.baatz.geomweight","Value for geometric (shape) weight");
        MandatoryOff("algorithm.baatz.geomweight");
        SetDefaultParameterFloat("algorithm.baatz.geomweight",  0.95);

        AddParameter(ParameterType_Int,"algorithm.meanshift.maxiter","max number of iterations");
        AddParameter(ParameterType_Float,"algorithm.meanshift.spatialr","Spatial bandwidth");
        MandatoryOff("algorithm.meanshift.spatialr");
        AddParameter(ParameterType_Float,"algorithm.meanshift.spectralr","Spectral bandwidth");
        MandatoryOff("algorithm.meanshift.spectralr");
        AddParameter(ParameterType_Float,"algorithm.meanshift.threshold","Threshold");
        MandatoryOff("algorithm.meanshift.threshold");
        AddParameter(ParameterType_Float,"algorithm.meanshift.ranger","Spectral range ramp");
        MandatoryOff("algorithm.meanshift.ranger");
        AddParameter(ParameterType_Choice,"algorithm.meanshift.modesearch","Activation of search mode");
        MandatoryOff("algorithm.meanshift.modesearch");
        AddChoice("algorithm.meanshift.modesearch.on","Activated");
        AddChoice("algorithm.meanshift.modesearch.off","Deactivated");

        // Processing Parameters
        AddParameter(ParameterType_Group,"processing","Set of parameters related to processing options");
        AddParameter(ParameterType_Int,"processing.memory","Maximum memory to be used on the main node");
        AddParameter(ParameterType_Int,"processing.nbproc","Number of processor allocated to the processing");
        MandatoryOff("processing.nbproc");
        AddParameter(ParameterType_Int,"processing.nbtilesperproc","Number tiles processed on each processor");
        MandatoryOff("processing.nbtilesperproc");
        AddParameter(ParameterType_Int,"processing.maxtilesizex","Maximum size of tiles along x axis");
        MandatoryOff("processing.maxtilesizex");
        AddParameter(ParameterType_Int,"processing.maxtilesizey","Maximum size of tiles along x axis");
        MandatoryOff("processing.maxtilesizey");
        AddParameter(ParameterType_Choice,"processing.writeimages","Activation of image traces");
        AddChoice("processing.writeimages.on","Activated");
        AddChoice("processing.writeimages.off","Deactivated");
        AddParameter(ParameterType_Choice,"processing.writegraphs","Activation of graph traces");
        AddChoice("processing.writegraphs.on","Activated");
        AddChoice("processing.writegraphs.off","Deactivated");
        AddParameter(ParameterType_Choice,"algorithm.baatz.aggregategraphs","Aggregation of graph traces");
        MandatoryOff("algorithm.baatz.aggregategraphs");
        AddChoice("algorithm.baatz.aggregategraphs.on","Activated");
        AddChoice("algorithm.baatz.aggregategraphs.off","Deactivated");

	// processing no data
	AddParameter(ParameterType_Float, "processing.nodatavalue", "Definition of no data value");
	SetDefaultParameterFloat("processing.nodatavalue",  0);	
	MandatoryOff("processing.nodatavalue");
    }

    void DoUpdateParameters()
    {
    	// If the user doesn't provide the tile size, retrieve if from the number of processor
		if(!HasUserValue("processing.maxtilesizex") || !HasUserValue("processing.maxtilesizey"))
		{
			if(HasUserValue("processing.nbproc") && HasUserValue("processing.nbtilesperproc") && HasUserValue("io.im") && HasUserValue("processing.memory"))
			{
				typedef otb::VectorImage<float, 2> ImageType;
				typedef otb::ImageFileReader<ImageType> ReaderType;
				ReaderType::Pointer reader = ReaderType::New();
				reader->SetFileName(GetParameterString("io.im"));
				reader->UpdateOutputInformation();
				uint32_t nb_band = reader->GetOutput()->GetNumberOfComponentsPerPixel();
				uint32_t width = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
				uint32_t height = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
				uint32_t nbproc = GetParameterInt("processing.nbproc");
				uint32_t tilesPerProc = GetParameterInt("processing.nbtilesperproc");
				uint32_t memPerProc = GetParameterInt("processing.memory");
				uint32_t maxIter = GetParameterInt("algorithm.baatz.maxiter");
				tilingConf = obia::SpatialTools::TilePartitionningOptimizer(nb_band, width, height, nbproc, tilesPerProc, memPerProc, maxIter);
				SetParameterInt("processing.maxtilesizex", tilingConf.tileWidth);
				SetParameterInt("processing.maxtilesizey", tilingConf.tileHeight);
			}
		}
		else if (tilingConf.tileWidth == 0 || tilingConf.tileHeight == 0)
		{
			tilingConf.tileWidth = GetParameterInt("processing.maxtilesizex");
			tilingConf.tileHeight = GetParameterInt("processing.maxtilesizey");
    	}
		if(HasUserValue("algorithm.baatz.numitfirstpartial") && tilingConf.maxIterPossible != 0)
		{
			if(GetParameterInt("algorithm.baatz.numitfirstpartial") > tilingConf.maxIterPossible)
			{
				otbAppLogWARNING("The number of iterations for first partial segmentation is higher than the maximum "
						"possible number of iteration according to the available memory. Setting it to " << tilingConf.maxIterPossible);
				SetParameterInt("algorithm.baatz.numitfirstpartial", tilingConf.maxIterPossible);
			}
		}
    }

    // Execute App
    void DoExecute()
    {
        
        /* Global parameters */
        std::string filename = GetParameterString("io.im");
        std::string outDir = GetParameterString("io.out.dir");
        std::string labelImage = GetParameterString("io.out.labelimage");
        std::string tmpDir = GetParameterString("io.temp");
        
        if(!HasUserValue("processing.maxtilesizex") || !HasUserValue("processing.maxtilesizey"))
        {
        	otbAppLogFATAL(<<"You need to provide either the tiles size, or the number of processor and the number of tiles per processor.");
        }
        uint32_t maxTileWidth = GetParameterInt("processing.maxtilesizex");
        uint32_t maxTileHeight = GetParameterInt("processing.maxtilesizey");
        unsigned long int memory = GetParameterInt("processing.memory");
        bool writeImages = false;
        bool writeGraphs = false;
        bool processUserNodata = false;

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

	if(HasUserValue("processing.nodatavalue"))
	{
		processUserNodata = true;
	}	        

	float noDataValue = GetParameterFloat("processing.nodatavalue");

        // BAATZ
        using InputImageType = otb::VectorImage<float, 2>;
        using LSBaatzSegmentationSchedulerType = otb::obia::LSBaatzSegmentationScheduler<InputImageType>;

        // MEANSHIFT
        using LabelPixelType              = unsigned int;
        using LSMeanShiftSchedulerType       = otb::obia::LSMeanShiftScheduler<InputImageType, LabelPixelType>;

        switch(GetParameterInt("algorithm"))
        {
            case ALG_BAATZ :
            {
                bool aggregateGraphs = false;
                switch (GetParameterInt("algorithm.baatz.aggregategraphs"))
                {
                    case ON:
                        aggregateGraphs = true;
                        break;
                    case OFF:
                        aggregateGraphs = false;
                        break;
                }
                uint32_t nbStartingIterations = GetParameterInt("algorithm.baatz.numitfirstpartial");
                uint32_t nbPartialIterations = GetParameterInt("algorithm.baatz.numitpartial");
		unsigned int maxIter = GetParameterInt("algorithm.baatz.maxiter");
                float threshold = GetParameterFloat("algorithm.baatz.scale");
                float decreasing = GetParameterFloat("algorithm.baatz.mindec");
                threshold = threshold * threshold;
                float spectralW = GetParameterFloat("algorithm.baatz.spectralweight");
                float shapeW = GetParameterFloat("algorithm.baatz.geomweight");

                auto lsBaatzFilter = LSBaatzSegmentationSchedulerType::New();
                lsBaatzFilter->SetFileName(filename);
                lsBaatzFilter->SetMaxTileSizeX(maxTileWidth);
                lsBaatzFilter->SetMaxTileSizeY(maxTileHeight);
                lsBaatzFilter->SetStartingNumberOfIterations(nbStartingIterations);
                lsBaatzFilter->SetPartialNumberOfIterations(nbPartialIterations);
		lsBaatzFilter->SetMaxNumberOfIterations(maxIter);
                lsBaatzFilter->SetTemporaryDirectory(tmpDir);
                lsBaatzFilter->SetAvailableMemory(memory);
                lsBaatzFilter->SetThreshold(threshold);
                lsBaatzFilter->SetDecreasing(decreasing);
                lsBaatzFilter->SetSpectralWeight(spectralW);
                lsBaatzFilter->SetShapeWeight(shapeW);
                lsBaatzFilter->SetWriteLabelImage(writeImages);
                lsBaatzFilter->SetWriteGraph(writeGraphs);
                lsBaatzFilter->SetAggregateGraphs(aggregateGraphs);
                lsBaatzFilter->SetOutputDir(outDir);
                lsBaatzFilter->SetLabelImageName(labelImage);
		lsBaatzFilter->SetProcessNoData(processUserNodata);
		lsBaatzFilter->SetNoDataValue(noDataValue);
	
                lsBaatzFilter->Update();

                break;
            }
            case ALG_MEANSHIFT : 
            {
                unsigned int maxIter = GetParameterInt("algorithm.meanshift.maxiter");
                unsigned int spatialr = GetParameterFloat("algorithm.meanshift.spatialr");
                float spectralr = GetParameterFloat("algorithm.meanshift.spectralr");
                float threshold = GetParameterFloat("algorithm.meanshift.threshold");
                float ranger = GetParameterFloat("algorithm.meanshift.ranger");
                                bool modeSearch = false;
                switch (GetParameterInt("algorithm.meanshift.modesearch")) {
                    case ON:
                        modeSearch = true;
                                        default:
                                                modeSearch = false;
                }
                                auto lsMSFilter = LSMeanShiftSchedulerType::New();
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
                lsMSFilter->SetWriteLabelImage(writeImages);
                lsMSFilter->SetWriteGraph(writeGraphs);
		lsMSFilter->SetProcessNoData(processUserNodata);
		lsMSFilter->SetNoDataValue(noDataValue);

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
