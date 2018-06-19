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
#ifndef otbObiaGraphToLabelImageFilter_txx
#define otbObiaGraphToLabelImageFilter_txx
#include "otbObiaGraphToLabelImageFilter.h"
#include "otbObiaConstExpr.h"

#include <unordered_set>

namespace otb
{
namespace obia
{



template< typename TInputGraph, typename TOutputImage>
void
GraphToLabelImageFilter<TInputGraph, TOutputImage>::
GenerateData()
{
    // Temporarily copy and label input graph
    auto graph = const_cast< InputGraphType * >( this->GetInput() );

    InputGraphPointerType labeledCopy = InputGraphType::New();
    labeledCopy->CopyGraph(graph);

    // Retrieve the dimension of the image to produce
    typename OutputImageType::IndexType index;
    typename OutputImageType::SizeType size;
    typename OutputImageType::RegionType region;
    typename OutputImageType::InternalPixelType noDataLabel = 0;

    auto image = this->GetOutput();
    image->Allocate();
    image->FillBuffer(noDataLabel);

    OutputImageIteratorType it(image, image->GetLargestPossibleRegion());
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        it.Set(0);
    }

    size[0] = labeledCopy->GetImageWidth();
    size[1] = labeledCopy->GetImageHeight();

    auto lambdaLabelWriter = [&](NodeType& node){
        auto label = node.m_Id + 1;
        std::unordered_set<CoordValueType> borderPixels;
        node.m_Contour.GenerateBorderPixels(borderPixels, size[0]);

        for(auto& pix : borderPixels){
            index[0] = pix % size[0];
            index[1] = pix / size[0];
            image->SetPixel(index, label);
        }
    };
    labeledCopy->ApplyForEachNode(lambdaLabelWriter);

    //Initialize reverse lut
    m_ReverseLut.clear();
    m_ReverseLut.resize(labeledCopy->GetNumberOfNodes()+1, noDataLabel);

    // Re-order node labels(left->right to top->bottom) (Remi Cresson correction)
    std::vector<typename OutputImageType::InternalPixelType> lut(labeledCopy->GetNumberOfNodes()+1, noDataLabel);
    typename OutputImageType::InternalPixelType label = 1;
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        auto inputLabel = it.Get();
        if(inputLabel != noDataLabel && lut[inputLabel] == noDataLabel)
          {
            lut[ inputLabel ] = label;
            m_ReverseLut[label] = inputLabel;
            label++;
          }
      }
    // Apply lut
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		it.Set(lut[it.Get()]);
	}
}


} // end of namespace obia
} // end of namespace otb

#endif
