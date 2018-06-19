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
        SetName("ComputeAttributes");
        SetDescription("Vector Attributes Computation Application");

        //Documentation
        SetDocName("Compute attributes");
        SetDocLongDescription("This application computes several attributes of a vector file");
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

    }

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

        //Compute attributes
        using ComputeAttributesFilterType = otb::obia::ComputeAttributesFilter<InputImageType>;
        auto computeAttributes = ComputeAttributesFilterType::New();

        const ogr::DataSource* outputDs =ogr::DataSource::New(vectorFile); //default 2nd argument is read
        computeAttributes->SetOGRData(outputDs);
        computeAttributes->SetInput(imgReader->GetOutput());
        computeAttributes->Update();
        std::cout << "End application" << std::endl;
    }
};

OTB_APPLICATION_EXPORT(GraphPolygonize)
}
}
