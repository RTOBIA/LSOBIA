#ifndef otbObiaSmallRegionsMergingFilter_txx
#define otbObiaSmallRegionsMergingFilter_txx
#include "otbObiaSmallRegionsMergingFilter.h"
namespace otb
{
namespace obia
{

template< typename TGraph >
SmallRegionsMergingFilter<TGraph>
::SmallRegionsMergingFilter() :
m_MinimalSurface(1), m_MergeOver(false), m_NumberOfIterations(-1)
{

}

template< typename TGraph >
SmallRegionsMergingFilter<TGraph>
::~SmallRegionsMergingFilter()
{
}

template< typename TGraph >
void
SmallRegionsMergingFilter<TGraph>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
	Superclass::PrintSelf(os, indent);
}

template< typename TGraph >
void
SmallRegionsMergingFilter<TGraph>
::GenerateData()
{
	std::cout << "SmallRegionsMergingFilter: Generate Data" << std::endl;
	std::cout << "Surface : " << m_MinimalSurface << std::endl;
	// Graft the input to the output
	// #TODO it may be good to create an abstract filter otbObiaInPlaceGraphFilter for this purpose
	// by following the same scheme of itkInPlaceImageFilter.
	GraphType * inputGraphPtr = dynamic_cast< GraphType * >( this->itk::ProcessObject::GetInput(0));

	auto outputGraph = this->GetOutput();
	outputGraph->GraftGraphByMove(inputGraphPtr);
	// Release the input graph.
	inputGraphPtr->Reset();

	//Iterate until no small regions found (meaning : m_MergeOver = true)
	m_MergeOver = false;
	uint32_t currIt = 0;
	if(m_NumberOfIterations < 0)
	{
		//Infinite loop
		while(!m_MergeOver)
		{
			std::cout << "Merge iteration : " << currIt << " size graph : " << outputGraph->GetNumberOfNodes() << std::endl;
			//IsGraphValid();
			m_MergeOver = DoOneIteration();
			currIt++;
		}
	}
	else
	{
		uint32_t current_iteration = 0;
		while(current_iteration < m_NumberOfIterations)
		{
			std::cout << "Iteration merge " << current_iteration << std::endl;
			current_iteration++;
			m_MergeOver = DoOneIteration();
		}
	}

	std::cout << "End Merging" << std::endl;
}

template< typename TGraph >
bool
SmallRegionsMergingFilter<TGraph>
::DoOneIteration()
{

	//Check minimal size : if size too small, look adjacent node and get the most similar in order to merge
	//Delete the merged node and update edge cost
	//Flag to indicate processing not over
	auto outputGraph = const_cast< GraphType * > (this->GetOutput() );
	bool merge_over = true;

	//Number of nodes
	uint32_t nbNodes = outputGraph->GetNumberOfNodes();
	//std::cout << "Parcours des noeuds : " << nbNodes << std::endl;

	//Looping across nodes
	//for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	for(uint32_t k = 0; k < nbNodes; ++k)
	{
		//std::cout << "Noeud " << k << std::endl;
		NodeType* nodeIt = outputGraph->GetNodeAt(k);
		if(!nodeIt->m_HasToBeRemoved)
		{
			//Check if node is too small
			if(nodeIt->m_Attributes.m_Area <= m_MinimalSurface)
			{

				//Get most similar adjacent node
				auto similar_node = GetMostSimilarNode(nodeIt);//GetMostSimilarNode(&(*nodeIt));

				//std::cout << "Small Area (" << nodeIt->m_Attributes.m_Area << " < " << m_MinimalSurface << ") for node " << nodeIt->m_Id << std::endl;
				//std::cout << "Target : " << similar_node->m_Id << std::endl;
				if(similar_node != nullptr)
				{
					auto nodeIn = nodeIt; //&(*nodeIt);
					auto nodeOut = similar_node;

					//Check ID, and eventually swap
					if(nodeIn->m_Contour.GetStartingCoords() > nodeOut->m_Contour.GetStartingCoords())
					{
						//std::cout << "SWAP " << nodeIn->m_Id << " and " << nodeOut->m_Id << std::endl;
						nodeIn = similar_node;
						nodeOut = nodeIt;//&(*nodeIt);
					}

					//Merge node
					//std::cout << "Merge NODE " << nodeIn->m_Id  << " with " << nodeOut->m_Id << std::endl;
					this->Merge(nodeIn, nodeOut);

					//Update attributes
					//std::cout << "Update attributes " << std::endl;
					UpdateSpecificAttributes(nodeIn, nodeOut);

					//Flag to indicate small region has been found (merge not over)
					merge_over = false;
				}

			}
		}
	}

	//Remove nodes who merged
	std::cout << "Remove Nodes" << std::endl;
	outputGraph->RemoveNodes();

	//Reconditionning graph
	auto lambdaValide = [](NodeType& node){ node.m_Valid = true;};
	outputGraph->ApplyForEachNode(lambdaValide);

	return merge_over;
}


template< typename TGraph >
typename SmallRegionsMergingFilter<TGraph>::NodeType*
SmallRegionsMergingFilter<TGraph>
::GetMostSimilarNode(NodeType * node)
{
	auto outputGraph = const_cast< GraphType * > (this->GetOutput() );
	double min_spec_distance = std::numeric_limits<double>::max();

	NodeType * similar_node = nullptr;

	//Looping across adjacent node
	//Get Edges
	for(auto edgeIt = node->Begin(); edgeIt != node->End(); edgeIt++)
	{
		//Get Target Node
		auto target_node = outputGraph->GetNodeAt(edgeIt->m_TargetId);

		if(target_node->m_Valid)
		{
			//Compute Spectral Distance
			double spec_distance(ComputeSpectralDistance(node, target_node));

			//If this distance is the minimum, update value
			if(spec_distance < min_spec_distance)
			{
				//Update distance
				min_spec_distance = spec_distance;

				//Update node
				similar_node = target_node;
			}
		}
	}

	//Look for the most similar node
	//Eventually swap the node and found node according to their id
	//Something like BaatzSegmentationFilter line 150
	return similar_node;

}

template< typename TGraph >
void
SmallRegionsMergingFilter<TGraph>
::UpdateMergingCosts()
{

}

template< typename TGraph >
void
SmallRegionsMergingFilter<TGraph>
::Merge(NodeType* nodeIn, NodeType * nodeOut)
{
	//UpdateSpecificAttributes(nodeIn, nodeOut);
	auto outputGraph = this->GetOutput();
	outputGraph->Merge(nodeIn, nodeOut);
}

template< typename TGraph >
void
SmallRegionsMergingFilter<TGraph>
::UpdateSpecificAttributes(NodeType * nodeIn, NodeType * nodeOut)
{
	//Update surface
	double area_1 = nodeIn->m_Attributes.m_Area;
	double area_2 = nodeOut->m_Attributes.m_Area;

	nodeIn->m_Attributes.m_Area = area_1 + area_2;
	//std::cout << "Old area " << area_1 << " New area : " << nodeIn->m_Attributes.m_Area << "(" << area_1 << " + " << area_2 << ")";

	//Update average spectral value
	//Number of bands
	uint32_t nbBands = nodeIn->m_Attributes.m_AvgSpec.size();

	for(uint32_t band = 0; band < nbBands; ++band)
	{
		float avgSpec_1 = nodeIn->m_Attributes.m_AvgSpec[band];
		float avgSpec_2 = nodeOut->m_Attributes.m_AvgSpec[band];

		nodeIn->m_Attributes.m_AvgSpec[band] = (area_1*avgSpec_1 + area_2*avgSpec_2)/(area_1 + area_2);
	}

	nodeIn->m_Valid = true;
}

template< typename TGraph >
float
SmallRegionsMergingFilter<TGraph>
::ComputeSpectralDistance(NodeType* node1, NodeType* node2)
{
	//Number of bands
	uint32_t nbBands = node1->m_Attributes.m_AvgSpec.size();

	//Distance
	double spec_distance = 0;

	for(uint32_t band = 0; band < nbBands; ++band)
	{
		//Average Spectral Value for each node
		float avgSpec_1 = node1->m_Attributes.m_AvgSpec[band];
		float avgSpec_2 = node2->m_Attributes.m_AvgSpec[band];

		//Compute distance
		spec_distance += (avgSpec_1 - avgSpec_2)*(avgSpec_1 - avgSpec_2);
	}

	return sqrt(spec_distance);
}

template< typename TGraph >
bool
SmallRegionsMergingFilter<TGraph>
::HasNoSmallRegions()
{
	auto outputGraph = const_cast< GraphType * > (this->GetOutput() );
	for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		if(nodeIt->m_Attributes.m_Area < this->m_MinimalSurface)
		{
			std::cout << "Graphe non valide!" << std::endl;
			exit(1);
		}
	}
}

template< typename TGraph >
bool
SmallRegionsMergingFilter<TGraph>
::IsGraphValid()
 {
	auto outputGraph = const_cast< GraphType * > (this->GetOutput() );
	for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	{
		if(!nodeIt->m_Valid)
		{
			std::cout << "Graphe non valide!" << std::endl;
			exit(1);
		}
	}
 }
}//End obia
}//End otb

#endif
