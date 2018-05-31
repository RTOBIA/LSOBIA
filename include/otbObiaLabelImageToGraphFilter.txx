#ifndef otbObiaLabelImageToGraphFilter_txx
#define otbObiaLabelImageToGraphFilter_txx
#include "otbObiaLabelImageToGraphFilter.h"

namespace otb
{
namespace obia
{

template< typename TLabelPixelType >
LabelImageToGraphFilter<TLabelPixelType>::
LabelImageToGraphFilter()
{
}

template< typename TLabelPixelType >
LabelImageToGraphFilter<TLabelPixelType>::
~LabelImageToGraphFilter()
{
}

template< typename TLabelPixelType >
void
LabelImageToGraphFilter<TLabelPixelType>::
GenerateData()
{
    // First step is to initialize the output graph
    InitOutput();

    // Second step is to build the graph
    BuildOutput();
}

template< typename TLabelPixelType >
void
LabelImageToGraphFilter<TLabelPixelType>::
BuildOutput()
{
    // It is an iterative process similar to a segmentation
    // process.
    while(DoOneIteration()){}
}

template< typename TLabelPixelType >
bool
LabelImageToGraphFilter<TLabelPixelType>::
DoOneIteration()
{
    // Get pointer to the output graph
    auto outputGraph = this->GetOutput();
    
    bool hasMerged = false;

    // Declare a lambda function
    auto lambdaDoIteration = [&](NodeType& node){

        // Determine if there is an adjacent node with the same label.
        auto adjNode = GetAdjacentNodeWithSameLabel(node);

        if(adjNode != nullptr)
        {
            outputGraph->Merge(&node, adjNode);

            node.m_Attributes.m_ListOfPixels.insert(node.m_Attributes.m_ListOfPixels.end(),
                                                    adjNode->m_Attributes.m_ListOfPixels.begin(),
                                                    adjNode->m_Attributes.m_ListOfPixels.end());

            hasMerged = true;
        } 
    };

    // Apply the lambda function to each node
    outputGraph->ApplyForEachNode(lambdaDoIteration);

    // Remove merged nodes
    outputGraph->RemoveNodes();

    return hasMerged;
}

template< typename TLabelPixelType >
typename LabelImageToGraphFilter<TLabelPixelType>::NodeType*
LabelImageToGraphFilter<TLabelPixelType>::GetAdjacentNodeWithSameLabel(NodeType& node)
{
    if(!node.m_HasToBeRemoved)
    {
        auto labelNode = node.m_Attributes.m_Label;
        auto outputGraph = this->GetOutput();

        // The lambda predicate
        auto predicate = [&labelNode, &outputGraph](const EdgeType& edg){

            auto adjNode = outputGraph->GetNodeAt(edg.m_TargetId);
            
            return (adjNode->m_Attributes.m_Label == labelNode && !adjNode->m_HasToBeRemoved);
        };

        auto edgeToAdjIt = node.FindEdgeIf(predicate);

        return ((edgeToAdjIt != node.m_Edges.end()) ? outputGraph->GetNodeAt(edgeToAdjIt->m_TargetId) : nullptr);
    }
    else
    {
        return nullptr;
    }
}

template< typename TLabelPixelType >
void
LabelImageToGraphFilter<TLabelPixelType>::
InitOutput()
{
    auto inputPtr = const_cast< InputImageType * >( this->GetInput() );
    auto outputGraph = this->GetOutput();
    
    // Retrieve the image dimensions
    const uint32_t imageWidth = inputPtr->GetLargestPossibleRegion().GetSize()[0];
    const uint32_t imageHeight = inputPtr->GetLargestPossibleRegion().GetSize()[1];

    outputGraph->SetImageWidth(imageWidth);
    outputGraph->SetImageHeight(imageHeight);
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
    	noDataFlags.reserve(1);
    	noDataValues.reserve(1);
    	noDataFlags.push_back(true);
    	noDataValues.push_back(value);
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

        // Specific initialization of the node and edge attributes
        // Here there are only attributes for the node
        newNode->m_Attributes.m_Label = it.Get();

        // Add the first pixel in the list of pixels
        newNode->m_Attributes.m_ListOfPixels.push_back(id);

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
                }
              }
            }

          }
        }
}


/**
 * Standard "PrintSelf" method
 */
template< typename TLabelPixelType >
void
LabelImageToGraphFilter<TLabelPixelType>::
PrintSelf(std::ostream & os, itk::Indent indent) const
{
     Superclass::PrintSelf(os, indent);
}

} // end of namespace obia
} // end of namespace otb

#endif
