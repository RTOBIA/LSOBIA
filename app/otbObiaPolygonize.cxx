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
#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbObiaGraph.h"
#include "otbObiaLabelImageToGraphFilter.h"
#include "otbObiaPolygonizeFilter.h"
#include "otbObiaImageToGraphFilter.h"
#include "otbObiaGraphOperations.h"
#include <string>
#include <sstream>

using std::string;
using std::stringstream;

namespace otb
{
namespace Wrapper
{

class Polygonize : public Application
{
public:

    typedef Polygonize Self;
    typedef Application SuperClass;
    typedef itk::SmartPointer<Self> Pointer;

    itkNewMacro(Self);
    itkTypeMacro(Polygonize, Application);

private:

    // Available input types
	enum InputType
	{
		GRAPH,
		LABEL_IMAGE
	};


    // Init App
    void DoInit()
    {

        //General description
        SetName("Polygonize");
        SetDescription("Polygonize Application");

        //Documentation
        SetDocName("Polygonization");
        SetDocLongDescription("This application provides several methods to perform polygonization of very high resolution images");
        SetDocLimitations("None");
        SetDocAuthors("OBIA-Team");
        SetDocSeeAlso(" ");

        // IO Parameters
        AddParameter(ParameterType_Group,"io","Set of parameters related to input/output");

			AddParameter(ParameterType_Choice,"io.input","The object to be polygonized");
			AddChoice("io.input.gr", "Input is a graph");
			AddChoice("io.input.img", "Input is a label image");

				AddParameter(ParameterType_String,  "io.input.gr.path", "Input graph path");
				SetParameterDescription("io.input.gr.path", "Path to the graph to be polygonized");

				AddParameter(ParameterType_String,  "io.input.img.path", "Input label image path");
				SetParameterDescription("io.input.img.path", "Path to the label image to be polygonized");

				AddParameter(ParameterType_Group, "io.out",  "Output directory");

					AddParameter(ParameterType_Directory, "io.out.dir",  "Output directory");
					SetParameterDescription("io.out.dir", "Output Directory");

					AddParameter(ParameterType_String, "io.out.gmlfile",  "GML FileName");
					SetParameterDescription("io.out.gmlfile", "GML FileName");
					MandatoryOff("io.out.gmlfile");

			AddParameter(ParameterType_Directory, "io.temp",  "Directory used for temporary data");
			SetParameterDescription("io.temp", "Temporary directory");

    }

    void DoUpdateParameters()
    {
    }

    // Execute App
    void DoExecute()
    {

        /* Global parameters */
        std::string outDir = GetParameterString("io.out.dir");
        std::string gmlFile = GetParameterString("io.out.gmlfile");
        std::string tmpDir = GetParameterString("io.temp");

        /* Useful definitions */
        using LabelPixelType = unsigned int;
        using LabelImageType = otb::Image< LabelPixelType, 2 >;
		using InputImageReader = otb::ImageFileReader< LabelImageType >;
		using LabelImageToGraphFilterType = otb::obia::LabelImageToGraphFilter< LabelPixelType >;
        using GraphType = typename LabelImageToGraphFilterType::OutputGraphType;

        /* Pipeline */

        // PolygonizeFilter
        using PolygonizeFilterType = otb::obia::PolygonizeFilter< GraphType >;
		auto polygonizeFilter = PolygonizeFilterType::New();

        switch(GetParameterInt("io.input"))
        {
			case GRAPH:
			{
				// Read the graph from disk
				auto graph = otb::obia::GraphOperations< GraphType >::ReadGraphFromDisk(GetParameterString("io.input.gr.path"));

				// Set the graph as input of the Polygonize filter
				polygonizeFilter->SetInput(graph);
				polygonizeFilter->Update();
				break;
			}
			case LABEL_IMAGE:
			{
				// Read the image

				InputImageReader::Pointer imageReader = InputImageReader::New();
				imageReader->SetFileName(GetParameterString("io.input.img.path"));

				// Convert the label image to graph
				LabelImageToGraphFilterType::Pointer labelImageToGraphFilter = LabelImageToGraphFilterType::New();
				labelImageToGraphFilter->SetInput(imageReader->GetOutput());
				labelImageToGraphFilter->Update();

				// Set the graph as input of the Polygonize filter
				polygonizeFilter->SetInput(labelImageToGraphFilter->GetOutput());
				break;
			}
			default:
			{
				otbAppLogFATAL(<< "No input provided");
				break;
			}
        }

        // Pipeline update
        polygonizeFilter->Update();
    }
};

OTB_APPLICATION_EXPORT(Polygonize)
}
}
