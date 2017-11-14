#ifndef otbObiaImageToBaatzGraphFilter_txx
#define otbObiaImageToBaatzGraphFilter_txx
#include "otbObiaImageToBaatzGraphFilter.h"
#include "otbObiaStreamUtils.h"

namespace otb
{
namespace obia
{

uint64_t
BaatzEdgeAttribute::GetMemorySize() const
{
    return FloatSize + BoolSize;
}

uint64_t 
BaatzEdgeAttribute::GetNumberOfBytesToSerialize() const
{
    // No specific data needs to be serialized
    return 0;
}

uint64_t 
BaatzNodeAttribute::GetMemorySize() const
{
    return m_AvgSpec.size() * 2 * FloatSize + // Mean and standard deviations
        2 * sizeof(std::vector<float>) +
      2 * UInt32Size + // area and perimeter
      BoolSize; // the flag
}

uint64_t 
BaatzNodeAttribute::GetNumberOfBytesToSerialize() const
{
  return stream_offset(m_AvgSpec) 
    + stream_offset(m_StdSpec) 
    + stream_offset(m_Area) 
    + stream_offset(m_Perimeter);
}

void 
BaatzNodeAttribute::Serialize(std::vector<char>& stream, uint64_t& position) const
{

    // Serialize the average spectral values
    to_stream(stream,m_AvgSpec,position);

    // Serialize the standard deviation values
    to_stream(stream,m_StdSpec,position);

    // Serialize the perimeter
    to_stream(stream,m_Perimeter,position);

    // Serialize the area
    to_stream(stream,m_Area,position);
}

void 
BaatzNodeAttribute::DeSerialize(const std::vector<char>& stream, uint64_t& position)
{
    // Deserialize the average spectral values
    from_stream(stream,m_AvgSpec,position);

    // Deserialize the standard deviation values
    from_stream(stream,m_StdSpec,position);

    // Deserialize the perimeter
    from_stream(stream,m_Perimeter,position);

    // Deserialize the area
    from_stream(stream,m_Area,position);
}

template< typename TInputImage >
ImageToBaatzGraphFilter< TInputImage >
::ImageToBaatzGraphFilter()
{
}

template< typename TInputImage >
ImageToBaatzGraphFilter< TInputImage >
::~ImageToBaatzGraphFilter()
{
}

template< typename TInputImage >
void
ImageToBaatzGraphFilter< TInputImage >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

template< typename TInputImage >
void
ImageToBaatzGraphFilter< TInputImage >
::GenerateData()
{
    auto inputPtr = const_cast< InputImageType * >( this->GetInput() );
    auto outputGraph = this->GetOutput();
    
    // Retrieve the image dimensions
    const uint32_t imageWidth = inputPtr->GetLargestPossibleRegion().GetSize()[0];
    const uint32_t imageHeight = inputPtr->GetLargestPossibleRegion().GetSize()[1];
    const uint32_t numBands = inputPtr->GetNumberOfComponentsPerPixel();

    outputGraph->SetImageWidth(imageWidth);
    outputGraph->SetImageHeight(imageHeight);
    outputGraph->SetNumberOfSpectralBands(numBands);
    outputGraph->SetProjectionRef(inputPtr->GetProjectionRef());

    //Set origin
    outputGraph->SetOriginX(inputPtr->GetOrigin()[0]);
    outputGraph->SetOriginY(inputPtr->GetOrigin()[1]);

    // Set the right number of starting nodes for memory allocation
    // 1 pixel = 1 node
    outputGraph->SetNumberOfNodes(imageWidth * imageHeight);

    // Create an image iterator
    ImageRegionConstIteratorType it(inputPtr, inputPtr->GetLargestPossibleRegion());

    // The running id
    IdType id = 0;

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        // Add a new node in the graph
        auto newNode = outputGraph->AddNode();

        // Generic intialization of the internal attributes of the node
        outputGraph->InitStartingNode(newNode, id);

        // Statistical attributes
        newNode->m_Attributes.m_StdSpec.assign(numBands, 0.0f);
        for(uint32_t b = 0; b < newNode->m_Attributes.m_StdSpec.size(); b++)
        {
            newNode->m_Attributes.m_AvgSpec.push_back(it.Get()[b]);
        }

        // Geometric attributes
        newNode->m_Attributes.m_Area = 1;
        newNode->m_Attributes.m_Perimeter = 4;

        // Force to compute all the costs at the first iteration
        newNode->m_Attributes.m_HasPreviouslyMerged = true;

        for(auto edgeIt = newNode->m_Edges.begin(); edgeIt != newNode->m_Edges.end(); edgeIt++)
        {
            edgeIt->m_Attributes.m_MergingCost = std::numeric_limits<float>::max();
            edgeIt->m_Attributes.m_CostUpdated = false;
        }

        // Increment the id
        id++;

    } // end for(it.GoToBegin(); !it.IsAtEnd(); ++it)
}

} // end of namespace obia
} // end of namespace otb

#endif
