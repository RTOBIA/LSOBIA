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
    std::cout<<"GenerateData de GraphToLabelImageFilter"<<std::endl;
    // Temporarily copy and label input graph
    auto graph = const_cast< InputGraphType * >( this->GetInput() );

    std::cout<<"    Copy graph"<<std::endl;
    InputGraphPointerType labeledCopy = InputGraphType::New();
    labeledCopy->CopyGraph(graph);

    // Retrieve the dimension of the image to produce
    typename OutputImageType::IndexType index;
    typename OutputImageType::SizeType size;
    typename OutputImageType::RegionType region;
    typename OutputImageType::InternalPixelType noDataLabel = 0;

    auto image = this->GetOutput();
    std::cout<<"    Allocate"<<std::endl;
    image->Allocate();
    std::cout<<"    Fill buffer"<<std::endl;
    image->FillBuffer(noDataLabel);

    std::cout<<"    Set 0"<<std::endl;
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
    std::cout<<"End apply for each node"<<std::endl;
    bool t = true;
    // Re-order node labels(left->right to top->bottom) (Remi Cresson correction)
    std::vector<typename OutputImageType::InternalPixelType> lut(labeledCopy->GetNumberOfNodes()+1, noDataLabel);
    typename OutputImageType::InternalPixelType label = 1;
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        if (t)
        {
            std::cout<<".";t=false;
        }
            auto inputLabel = it.Get();
        if(inputLabel != noDataLabel && lut[inputLabel] == noDataLabel)
          {
            lut[ inputLabel ] = label;
            label++;
          }
      }
    std::cout<<std::endl;
    // Apply lut
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        it.Set(lut[it.Get()]);
      }
}


} // end of namespace obia
} // end of namespace otb

#endif
