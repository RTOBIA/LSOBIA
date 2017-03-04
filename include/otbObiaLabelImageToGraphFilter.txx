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

		// Increment the id
		id++;

	} // end for(it.GoToBegin(); !it.IsAtEnd(); ++it)
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
