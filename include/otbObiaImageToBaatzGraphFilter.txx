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
  return sizeof(BaatzEdgeAttribute);
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
  return sizeof(BaatzNodeAttribute)+ m_AvgSpec.size() * 2 * FloatSize;
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
  std::cout << "Generate Baatz graph" << std::endl;

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

  // Get no-data values
  bool noDataPresent = false;
  std::vector<bool> noDataFlags;
  std::vector<double> noDataValues;

  // default : no data indicated in metadata
  if (!this->m_ProcessNoData)
  {
	
  	noDataPresent = otb::ReadNoDataFlags(inputPtr->GetMetaDataDictionary(), noDataFlags, noDataValues);
  }
  else // no data defined by user
  {
	noDataPresent = true;
	double value = static_cast<double>(this->m_NoDataValue);
	noDataFlags.reserve(numBands);
	noDataValues.reserve(numBands);

	for (int i=0; i<numBands; i++)
	{
		noDataFlags.push_back(true);
		noDataValues.push_back(value);
	}	
  }

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

    // If the pixel is no-data, then the node has to be removed
    newNode->m_HasToBeRemoved = (noDataPresent && otb::IsNoData<double>(it.Get(), noDataFlags, noDataValues));

    // Increment the id
    id++;

    } // end for(it.GoToBegin(); !it.IsAtEnd(); ++it)


  // Add the edges
  for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {
    if (!nodeIt->m_HasToBeRemoved)
      {
      // Count the required number of edges
      auto neighbors = otb::obia::SpatialTools::FourConnectivity(nodeIt->m_Id, imageWidth, imageHeight);
      uint32_t numEdges = 0;
      for(unsigned short n = 0; n < 4; n++)
        {
        if(neighbors[n] > -1)
          if  (!outputGraph->GetNodeAt(neighbors[n])->m_HasToBeRemoved)
            numEdges++;

        }

      // Initialisation of the edges
      nodeIt->m_Edges.reserve(numEdges);
      for(unsigned short n = 0; n < 4; n++)
        {
        if(neighbors[n] > -1)
          {
          if  (!outputGraph->GetNodeAt(neighbors[n])->m_HasToBeRemoved)
            {
            // Add an edge to the current node targeting the adjacent node
            auto newEdge = nodeIt->AddEdge();

            // Add the target
            newEdge->m_TargetId = neighbors[n];

            // Initialisation of the boundary
            newEdge->m_Boundary = 1;

            // Initialisation of the merging cost
            newEdge->m_Attributes.m_MergingCost = std::numeric_limits<float>::max();
            }
          }
        }

      }
    }

  // Here, we remove the "no-data" nodes
  outputGraph->RemoveNodes();

  std::cout << "Nombre noeuds : " << outputGraph->GetNumberOfNodes() << std::endl;
  std::cout << "Adresse : " << outputGraph << std::endl;
}

} // end of namespace obia
} // end of namespace otb

#endif
