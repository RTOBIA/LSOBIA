#ifndef otbObiaImageToBaatzGraphFilter_txx
#define otbObiaImageToBaatzGraphFilter_txx
#include "otbObiaImageToBaatzGraphFilter.h"

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
    // The flag does not need to be serialized.
    return UInt32Size + 
    m_AvgSpec.size() * 2 * FloatSize +
    2* UInt32Size;
}

void 
BaatzNodeAttribute::Serialize(std::vector<char>& stream, uint64_t& position) const
{
    // Serialize the number of bands
    const uint32_t numBands = m_AvgSpec.size();
    std::memcpy(&(stream[position]), &numBands, UInt32Size);
    position += UInt32Size;

    // Serialize the average spectral values
    std::memcpy(&(stream[position]), &m_AvgSpec[0], numBands * FloatSize);
    position += numBands * FloatSize;

    // Serialize the standard deviation values
    std::memcpy(&(stream[position]), &m_StdSpec[0], numBands * FloatSize);
    position += numBands * FloatSize;

    // Serialize the perimeter
    std::memcpy(&(stream[position]), &m_Perimeter, UInt32Size);
    position += UInt32Size;

    // Serialize the area
    std::memcpy(&(stream[position]), &m_Area, UInt32Size);
    position += UInt32Size;
}

void 
BaatzNodeAttribute::DeSerialize(const std::vector<char>& stream, uint64_t& position)
{
    // Deserialize the number of bands.
    uint32_t numBands = 0;
    std::memcpy(&numBands, &(stream[position]), UInt32Size);
    position += UInt32Size;

    // Initialization of the statistical attributes
    m_AvgSpec.assign(numBands, 0.0f);
    m_StdSpec.assign(numBands, 0.0f);

    // Deserialize the average spectral values
    std::memcpy(&m_AvgSpec[0], &(stream[position]), numBands * FloatSize);
    position += numBands * FloatSize;

    // Deserialize the standard deviation values
    std::memcpy(&m_StdSpec[0], &(stream[position]), numBands * FloatSize);
    position += numBands * FloatSize;

    // Deserialize the perimeter
    std::memcpy(&m_Perimeter, &(stream[position]), UInt32Size);
    position += UInt32Size;

    // Deserialize the area
    std::memcpy(&m_Area, &(stream[position]), UInt32Size);
    position += UInt32Size;
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
  std::cout << "Generate Baatz graphe" << std::endl;
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
  std::vector<bool> noDataFlags;
  std::vector<double> noDataValues;
  bool noDataPresent = otb::ReadNoDataFlags(inputPtr->GetMetaDataDictionary(), noDataFlags, noDataValues);

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
