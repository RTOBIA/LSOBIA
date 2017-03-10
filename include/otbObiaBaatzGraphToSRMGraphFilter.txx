#ifndef otbObiaBaatzGraphToSRMGraphFilter_txx
#define otbObiaBaatzGraphToSRMGraphFilter_txx
#include "otbObiaBaatzGraphToSRMGraphFilter.h"

namespace otb
{
namespace obia
{

BaatzToSRMGraphFilter
::BaatzToSRMGraphFilter()
{
	std::cout << "BaatzToSRMGraphFilter: Convert" << std::endl;
	// Modify superclass default values, can be overridden by subclasses
  	this->SetNumberOfRequiredInputs(1);
}

BaatzToSRMGraphFilter
::~BaatzToSRMGraphFilter()
{
}

void
BaatzToSRMGraphFilter
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

void
BaatzToSRMGraphFilter
::GenerateData()
{
	std::cout << "BaatzToSRMGraphFilter: Generate Data" << std::endl;
	//The Baatz Graph already contains the mean and the surface of each node
	//No need a lot of work to convert this graph
	auto inputPtr 	 = const_cast< InputGraphType * >( this->GetInput() );
	auto outputGraph = const_cast< OutputGraphType * >(this->GetOutput() );

	// Retrieve the image dimensions
	outputGraph->SetImageWidth(inputPtr->GetImageWidth());
	outputGraph->SetImageHeight(inputPtr->GetImageHeight());
	outputGraph->SetNumberOfSpectralBands(inputPtr->GetNumberOfSpectralBands());
	outputGraph->SetProjectionRef(inputPtr->GetProjectionRef());

	// Set the right number of starting nodes for memory allocation
	// 1 pixel = 1 node
	outputGraph->SetNumberOfNodes(inputPtr->GetNumberOfNodes());

	//Looping across the nodes
	for(auto nodeIt = inputPtr->Begin(); nodeIt != inputPtr->End(); nodeIt++)
	{
		// Create a new node
		auto new_node = convertInputNode(&(*nodeIt));
		// Nothing is done on the output
		(void) new_node;
	}

	//Release input graph
	inputPtr->Reset();
}

typename BaatzToSRMGraphFilter::OutputNodeType*
BaatzToSRMGraphFilter
::convertInputNode(InputNodeType* node)
{
	auto outputGraph = const_cast< OutputGraphType * >(this->GetOutput() );
	OutputNodeType* output_node = outputGraph->AddNode();

	//Copy node attributes
	// Flag indicating if the node has to be removed from the graph
	output_node->m_HasToBeRemoved 		= node->m_HasToBeRemoved;
	output_node->m_Valid 		  		= node->m_Valid;
	output_node->m_Id			  		= node->m_Id;
	output_node->m_BoundingBox 	  		= node->m_BoundingBox;
	output_node->m_Contour		 	 	= node->m_Contour; //Check if copy or not
	output_node->m_Attributes.m_Area 	= node->m_Attributes.m_Area;
	output_node->m_Attributes.m_AvgSpec = node->m_Attributes.m_AvgSpec;
	output_node->m_Attributes.m_StdSpec = node->m_Attributes.m_StdSpec;

	//Get Edges
	for(auto edgeIt = node->Begin(); edgeIt != node->End(); edgeIt++)
	{
		OutputNodeType::EdgeType edge;

		// The position of the edge in the adjacency graph
		edge.m_TargetId = edgeIt->m_TargetId;

		// The boundary length between both adjacent nodes
		edge.m_Boundary = edgeIt->m_Boundary;

		//Add in the vector
		output_node->m_Edges.push_back(edge);
	}

	return output_node;

}
}
}

#endif
