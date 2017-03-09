#ifndef __otbObiaGenericRegionMergingFilter_txx
#define __otbObiaGenericRegionMergingFilter_txx
#include "otbObiaGenericRegionMergingFilter.h"

namespace otb
{
namespace obia
{

template< typename TInputGraph, 
		  typename TOutputGraph, 
		  typename TMergingCostFunc, 
		  typename THeuristic, 
		  typename TUpdateAttributeFunc>
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>
::GenericRegionMergingFilter()
: m_MaxNumberOfIterations(75), m_AppliedNumberOfIterations(0), m_MergingOver(false), m_MergingCostFunc(nullptr),
m_HeuristicFunc(nullptr), m_UpdateAttributeFunc(nullptr)
{
	std::cout << "Create Filter Object" << std::endl;
}

template< typename TInputGraph, 
		  typename TOutputGraph, 
		  typename TMergingCostFunc, 
		  typename THeuristic, 
		  typename TUpdateAttributeFunc>
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>
::~GenericRegionMergingFilter()
{
	if(m_MergingCostFunc != nullptr){
		delete m_MergingCostFunc;
	}

	if(m_HeuristicFunc != nullptr){
		delete m_HeuristicFunc;
	}

	if(m_UpdateAttributeFunc != nullptr){
		delete m_UpdateAttributeFunc;
	}
}


template< typename TInputGraph,
		  typename TOutputGraph,
		  typename TMergingCostFunc,
		  typename THeuristic,
		  typename TUpdateAttributeFunc>
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template< typename TInputGraph,
		  typename TOutputGraph,
		  typename TMergingCostFunc,
		  typename THeuristic,
		  typename TUpdateAttributeFunc>
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::
GenerateData()
{
	std::cout << "Generate Data" << std::endl;
	auto outputGraph = this->GetOutputByMove();

	//Set the graph for the heuristic
	this->GetHeuristicFunc()->SetGraph(outputGraph);

	// Reconditionning of the graph
	for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		nodeIt->m_HasToBeRemoved = false;
		nodeIt->m_Valid = true;
		nodeIt->m_Attributes.m_HasPreviouslyMerged = true;
		for(auto& edg : nodeIt->m_Edges)
		{
			edg.m_Attributes.m_MergingCost = m_MergingCostFunc->GetMax();
			edg.m_Attributes.m_CostUpdated = false;
		}
	}

	for(uint32_t i = 0; i < m_MaxNumberOfIterations; i++)
	{
		std::cout << "iteration " << i+1 << "/" << m_MaxNumberOfIterations <<  std::endl;
		if(!DoOneIteration())
		{
			m_MergingOver = true;
			break;
		}
	}
}


template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc>
bool 
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::DoOneIteration()
{
	// Compute the merging costs for all the pairs of adjacent nodes.
	ComputeMergingCosts();

	auto outputGraph = this->GetOutput();
	
	// Flag indicating if at least one merge has been done during the iteration.
	bool merged = false;

	for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		// Heuristic to determine with which adjacent node this current node has to merge.
		auto nodeIn = m_HeuristicFunc->GetBestAdjacentNode(&(*nodeIt));

		// The heuristic must return true if no adjacent node has been found.
		if(nodeIn != nullptr)
		{
			auto nodeOut = outputGraph->GetNodeAt(nodeIn->m_Edges.front().m_TargetId);

			auto cost = nodeIn->m_Edges.front().m_Attributes.m_MergingCost;
			
			Merge(nodeIn, nodeOut);

			merged = true;

		}

	}

	outputGraph->RemoveNodes();

	if(outputGraph->GetNumberOfNodes() < 2)
	{
		return false;
	}

	return merged;

}

template< typename TInputGraph, 
		  typename TOutputGraph, 
		  typename TMergingCostFunc, 
		  typename THeuristic,
		  typename TUpdateAttributeFunc>
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>
::Merge(NodeType* nodeIn, NodeType * nodeOut)
{
	m_UpdateAttributeFunc->UpdateAttributes(nodeIn, nodeOut);
	auto outputGraph = this->GetOutput();
	outputGraph->Merge(nodeIn, nodeOut);
}

template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc >
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::ComputeMergingCosts()
{
	std::cout << "Compute Merging cost" << std::endl;
	// Retrieve the output graph.
	auto outputGraph = this->GetOutput();


	MergingCostValueType minCost;
	uint64_t minNodeId, idx, minIdx;

	// The nodes must have a boolean attribute called m_CostUpdated.
	for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		for(auto edgeIt = nodeIt->m_Edges.begin(); edgeIt != nodeIt->m_Edges.end(); edgeIt++)
		{
			edgeIt->m_Attributes.m_CostUpdated = false;
		}
	}

	// Loop over the nodes
	for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		if(m_MergingCostFunc->ComputeMergingCostsForThisNode(&(*nodeIt)))
		{

			// The merging cost function must give the maximum value of the merging cost.
			minCost = MergingCostFunctionType::Max();
			idx = 0;
			minIdx = 0;
			nodeIt->m_HasToBeRemoved = false;
			nodeIt->m_Valid = true;

			// Loop over the edges
			for(auto edgeIt = nodeIt->m_Edges.begin(); edgeIt != nodeIt->m_Edges.end(); edgeIt++)
			{

				// Retrieve the adjacent node.
				auto adjNode = outputGraph->GetNodeAt(edgeIt->m_TargetId);

				if(m_MergingCostFunc->ComputeMergingCostsForThisAdjNode(adjNode))
				{

					// If the cost is not updated and if one of the adjacent nodes
					// has merged at the previous iteration then we must compute the
					// merging cost.
					if(!(edgeIt->m_Attributes.m_CostUpdated) && 
						(nodeIt->m_Attributes.m_HasPreviouslyMerged || 
					   	adjNode->m_Attributes.m_HasPreviouslyMerged))
					{
						edgeIt->m_Attributes.m_MergingCost = m_MergingCostFunc->ComputeMergingCost(&(*nodeIt), adjNode);
						auto adjNodeToCurr = adjNode->FindEdge(nodeIt->m_Id);
						adjNodeToCurr->m_Attributes.m_MergingCost = edgeIt->m_Attributes.m_MergingCost;
						edgeIt->m_Attributes.m_CostUpdated = true;
						adjNodeToCurr->m_Attributes.m_CostUpdated = true;
					}


					// If the current cost is minimum than we record it.
					if(edgeIt->m_Attributes.m_MergingCost < minCost)
					{
						minCost = edgeIt->m_Attributes.m_MergingCost;
						minNodeId = adjNode->GetFirstPixelCoords();
						minIdx = idx;
					}

					// In case of equality, we keep the adjacent node with the lower starting
					// coordinates.
					else if(minCost == edgeIt->m_Attributes.m_MergingCost)
					{
						if(adjNode->GetFirstPixelCoords() < minNodeId)
						{
							minNodeId = adjNode->GetFirstPixelCoords();
							minIdx = idx;
						}
					}

				} // end if(MergingCostFunctionType::ComputeMergingCostsForThisAdjNode(adjNode))

				idx++;

			} // end loop over the edges.

			// Finally we move the adjacent node with the lower merging cost
			// at the first position in the list of adjacent nodes.
			std::swap(nodeIt->m_Edges[0], nodeIt->m_Edges[minIdx]);

		} // end if(MergingCostFunctionType::ComputeMergingCostsForThisNode(&(*nodeIt)))

	} // end for loop over the nodes

	// Reset the merge flag for all the nodes.
	for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		nodeIt->m_Attributes.m_HasPreviouslyMerged = false;
	}

}

template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc >
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::CheckValidity()
{
	if(m_MergingCostFunc == nullptr || m_HeuristicFunc == nullptr || m_UpdateAttributeFunc == nullptr)
	{
		std::cerr << "GenericRegionMergingFilter not initialized like it should..." << std::endl;
	}
}
} // end of namespace obia
} // end of namespace otb

#endif 
