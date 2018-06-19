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
#ifndef otbObiaGraphToLabelImageFilter_h
#define otbObiaGraphToLabelImageFilter_h
#include "itkImageRegionIterator.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "otbObiaGraphToImageFilter.h"


namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputImage>
class GraphToLabelImageFilter : public GraphToImageFilter<TInputGraph, TOutputImage>
{
public:

    /** Some convenient alias */
    using Self                       = GraphToLabelImageFilter;
    using Superclass                 = itk::ImageSource<TOutputImage>;
    using Pointer                    = itk::SmartPointer<Self>;
    using ConstPointer               = itk::SmartPointer< const Self >;
    using InputGraphType             = TInputGraph;
    using InputGraphPointerType  	 = typename TInputGraph::Pointer;
    using NodeType                	 = typename InputGraphType::NodeType;
    using OutputImageType         	 = TOutputImage;
    using OutputImagePointerType  	 = typename TOutputImage::Pointer;
    using OutputImageIteratorType  	 = itk::ImageRegionIterator<OutputImageType>;
    using FillholeFilterType      	 = itk::GrayscaleFillholeImageFilter<OutputImageType,OutputImageType>;

    itkNewMacro(Self);

    //Get LUT
    std::vector<typename OutputImageType::InternalPixelType> GetReverseLut(){return m_ReverseLut;};

protected:

    virtual void GenerateData();

    //Keep LUT to get correspondance between label and node id
    std::vector<typename OutputImageType::InternalPixelType> m_ReverseLut;
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraphToLabelImageFilter.txx"
#endif
