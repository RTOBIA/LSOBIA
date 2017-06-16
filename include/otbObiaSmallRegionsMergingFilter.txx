#ifndef otbObiaSmallRegionsMergingFilter_txx
#define otbObiaSmallRegionsMergingFilter_txx

#include "otbObiaSmallRegionsMergingFilter.h"

namespace otb
{
namespace obia
{


template< typename TCost, typename TGraph >
SRMMergingCost<TCost, TGraph>
::SRMMergingCost() : m_MinimalSurface(0)
{
    std::cout << "New SMR Cost" << std::endl;
}

template< typename TCost, typename TGraph >
SRMMergingCost<TCost, TGraph>
::~SRMMergingCost()
{
}

template< typename TCost, typename TGraph >
bool
SRMMergingCost<TCost, TGraph>::ComputeMergingCostsForThisNode(NodeType* curNode)
{
	return true;
}

template< typename TCost, typename TGraph >
bool
SRMMergingCost<TCost, TGraph>::ComputeMergingCostsForThisAdjNode(NodeType* curNode)
{
	return true;
}

template< typename TCost, typename TGraph >
typename SRMMergingCost< TCost, TGraph>::ValueType
SRMMergingCost<TCost, TGraph>::ComputeMergingCost(NodeType* n1, NodeType* n2)
{
	if(n1->m_Attributes.m_Area > m_MinimalSurface &&
	   n2->m_Attributes.m_Area > m_MinimalSurface)
	{
		return std::numeric_limits<TCost>::max();//Return max distance
	}

	//Number of bands
	uint32_t nbBands = n1->m_Attributes.m_AvgSpec.size();

    //Distance
    double spec_distance = 0;

    for(uint32_t band = 0; band < nbBands; ++band)
    {
        //Average Spectral Value for each node
        float avgSpec_1 = n1->m_Attributes.m_AvgSpec[band];
        float avgSpec_2 = n2->m_Attributes.m_AvgSpec[band];

        //Compute distance
        spec_distance += (avgSpec_1 - avgSpec_2)*(avgSpec_1 - avgSpec_2);
    }

    return spec_distance;
}

template< typename TGraph >
SRMHeuristic<TGraph>
::SRMHeuristic()
{
}

template< typename TGraph >
SRMHeuristic<TGraph>
::~SRMHeuristic()
{
}
template< typename TGraph >
typename SRMHeuristic<TGraph>::NodeType*
SRMHeuristic<TGraph>
::GetBestAdjacentNode(NodeType* node)
{

	//If this node has an area too big, do not consider it
	if(node->m_Attributes.m_Area > m_MinimalSurface || !node->m_Valid){
		//std::cout << "Noeud " << node->m_Id << " a une surface de " << node->m_Attributes.m_Area << std::endl;
		return nullptr;
	}
	else
	{
		return node;
	}

	auto outputGraph = this->m_Graph;

	if(node->m_Valid)
	{
		//Get the first edge (lowest cost)
		//auto bestAdjNode = outputGraph->GetNodeAt(node->m_Edges.front().m_TargetId);
		auto bestAdjNode = GetMostSimilarNode(node);

		if(bestAdjNode != nullptr)
		{
			if(bestAdjNode->m_Valid)
			{
				if(CheckMutuality(node, bestAdjNode))
				{
					return bestAdjNode;

				}
			}
		}
	}

	return nullptr;

}

template< typename TGraph >
typename SRMHeuristic<TGraph>::NodeType*
SRMHeuristic<TGraph>
::GetBestSmallAdjacentNode(NodeType* nodeIn)
{

	auto outputGraph = this->m_Graph;
	NodeType* bestNode = nullptr;
	float minDistance = std::numeric_limits<float>::max();

	/**Go through the edge*/
	//Loop over edge
	for(auto edgeIt = nodeIn->Begin(); edgeIt != nodeIn->End(); edgeIt++)
	{
		auto adjNode = outputGraph->GetNodeAt(edgeIt->m_TargetId);
		if(adjNode->m_Attributes.m_Area <= m_MinimalSurface)
		{
			//std::cout << "Noeud adjacent a une surface " << adjNode->m_Attributes.m_Area << std::endl;
			//Compute distance
			float curDistance = ComputeSpectralDistance(nodeIn, adjNode);

			//Update minDistance
			if(curDistance < minDistance)
			{
				//std::cout << "Update distance" << std::endl;
				//Update min distance
				minDistance = curDistance;
				bestNode = adjNode;
			}else if(minDistance == curDistance)
			{
				//Get the node with the lowest ID
				if(adjNode->GetFirstPixelCoords() < bestNode->GetFirstPixelCoords())
				{
					//Update bestNode
					bestNode = adjNode;
				}
			}
		}
	}

	if(!bestNode->m_Valid)
	{
		return nullptr;
	}

	return bestNode;
}


template< typename TGraph >
typename SRMHeuristic<TGraph>::NodeType*
SRMHeuristic<TGraph>
::GetMostSimilarNode(NodeType* node)
 {
	auto outputGraph = this->m_Graph;
	NodeType* bestNode = nullptr;
	float minDistance = std::numeric_limits<float>::max();

	/**Go through the edge*/
	//Loop over edge
	for(auto edgeIt = node->Begin(); edgeIt != node->End(); edgeIt++)
	{
		auto adjNode = outputGraph->GetNodeAt(edgeIt->m_TargetId);
		//std::cout << "Noeud adjacent a une surface " << adjNode->m_Attributes.m_Area << std::endl;
		//Compute distance
		float curDistance = ComputeSpectralDistance(node, adjNode);

		//Update minDistance
		if(curDistance < minDistance)
		{
			//std::cout << "Update distance" << std::endl;
			//Update min distance
			minDistance = curDistance;
			bestNode = adjNode;
		}else if(minDistance == curDistance)
		{
			//Get the node with the lowest ID
			if(adjNode->GetFirstPixelCoords() < bestNode->GetFirstPixelCoords())
			{
				//Update bestNode
				bestNode = adjNode;
			}
		}
	}

	return bestNode;
 }

template< typename TGraph >
bool
SRMHeuristic<TGraph>
::CheckMutuality(NodeType* node_1, NodeType* node_2)
{
	auto outputGraph = this->m_Graph;
	float minDistance = std::numeric_limits<float>::max();
	NodeType* bestNode = nullptr;

	//Check if node_2 got node_1 for best adjacent node
	//We assume that node_1 got a surface < minimal_surface
	bestNode = GetBestSmallAdjacentNode(node_2);
	if(bestNode == nullptr)
	{
		return false;
	}

	if(bestNode->m_Id == node_1->m_Id)
	{
		return true;
	}
	else
	{
		return false;
	}

}

template< typename TGraph >
float
SRMHeuristic<TGraph>
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

    return spec_distance;
}
template< typename TGraph >
SRMUpdateAttribute<TGraph>
::SRMUpdateAttribute()
{

}

template< typename TGraph >
SRMUpdateAttribute<TGraph>
::~SRMUpdateAttribute()
{

}
template<typename TGraph >
void
SRMUpdateAttribute<TGraph>
::UpdateAttributes(NodeType * nodeIn, NodeType *  nodeOut)
 {

    //nodeIn->m_Valid = true;

    // Update statistical attributes
        const float aArea = static_cast<float>(nodeIn->m_Attributes.m_Area);
        const float bArea = static_cast<float>(nodeOut->m_Attributes.m_Area);
        const float a_sum = aArea + bArea;

        float mean;
        float stddev;
        float stddevTmp;

        for(unsigned int b = 0; b < nodeIn->m_Attributes.m_AvgSpec.size(); ++b)
        {
            mean = (nodeIn->m_Attributes.m_AvgSpec[b] * aArea + nodeOut->m_Attributes.m_AvgSpec[b] * bArea) / a_sum;
            stddevTmp = aArea * nodeIn->m_Attributes.m_StdSpec[b] * nodeIn->m_Attributes.m_StdSpec[b];
            stddevTmp += bArea * nodeOut->m_Attributes.m_StdSpec[b] * nodeOut->m_Attributes.m_StdSpec[b];
            stddevTmp += aArea * (nodeIn->m_Attributes.m_AvgSpec[b] - mean) * (nodeIn->m_Attributes.m_AvgSpec[b] - mean);
            stddevTmp += bArea * (nodeOut->m_Attributes.m_AvgSpec[b] - mean) * (nodeOut->m_Attributes.m_AvgSpec[b] - mean);
            stddev = std::sqrt(stddevTmp / a_sum);
            nodeIn->m_Attributes.m_AvgSpec[b] = mean;
            nodeIn->m_Attributes.m_StdSpec[b] = stddev;
        }

        nodeIn->m_Attributes.m_Area = a_sum;
        nodeIn->m_Attributes.m_Perimeter += (nodeOut->m_Attributes.m_Perimeter - 2*(nodeIn->FindEdge(nodeOut->m_Id)->m_Boundary));

        nodeIn->m_Attributes.m_HasPreviouslyMerged = true;
}


template< typename TGraph >
SmallRegionsMergingFilter<TGraph>
::SmallRegionsMergingFilter() :
m_MinimalSurface(1), m_MergingOver(false), m_NumberOfIterations(-1)
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

    //Iterate until no small regions found (meaning : m_MergingOver = true)
    m_MergingOver = false;
    uint32_t currIt = 0;
    if(m_NumberOfIterations < 0)
    {
        //Infinite loop
        while(!m_MergingOver)
        {
            std::cout << "Merge iteration : " << currIt << " size graph : " << outputGraph->GetNumberOfNodes() << std::endl;
            //IsGraphValid();
            m_MergingOver = DoOneIteration();
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
            m_MergingOver = DoOneIteration();
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
	for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
	//for(uint32_t k = 0; k < nbNodes; ++k)
	{
		auto node = nodeIt;
		//std::cout << "Noeud " << k << std::endl;
		//NodeType* nodeIt = outputGraph->GetNodeAt(k);

		if(!nodeIt->m_HasToBeRemoved)
		{
			//Check if node is too small
			if(nodeIt->m_Attributes.m_Area <= m_MinimalSurface)
			{

				//Get most similar adjacent node
				auto similar_node = GetMostSimilarNode(&(*node));//GetMostSimilarNode(&(*nodeIt));

				if(similar_node != nullptr && similar_node->m_Valid)
				{
					//Check mutuality
					if(CheckMutuality(&(*node), similar_node))
					{
						auto nodeIn = &(*node); //&(*nodeIt);
						auto nodeOut = similar_node;

						//Check ID, and eventually swap
						if(nodeIn->GetFirstPixelCoords() > nodeOut->GetFirstPixelCoords())
						{
							auto tmp = nodeOut;
							nodeOut = nodeIn;
							nodeIn = tmp;//&(*nodeIt);
						}

						//std::cout << "Merge NODE " << nodeIn->m_Id  << " with " << nodeOut->m_Id << std::endl;
						this->Merge(nodeIn, nodeOut);

						//Flag to indicate small region has been found (merge not over)
						merge_over = false;
					}
				}
			}
		}
	}

	//Remove nodes who merged
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
        //if(!target_node->m_HasToBeRemoved)
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
    UpdateSpecificAttributes(nodeIn, nodeOut);
    auto outputGraph = this->GetOutput();
    outputGraph->Merge(nodeIn, nodeOut);

    //nodeIn->m_Valid = true;

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
    /*std::cout << "Old area " << area_1 << " New area : "
              << nodeIn->m_Attributes.m_Area << "(" << area_1 << " + " << area_2 << ")" << std::endl;*/

    //Update average spectral value
    //Number of bands
    uint32_t nbBands = nodeIn->m_Attributes.m_AvgSpec.size();

    for(uint32_t band = 0; band < nbBands; ++band)
    {
        float avgSpec_1 = nodeIn->m_Attributes.m_AvgSpec[band];
        float avgSpec_2 = nodeOut->m_Attributes.m_AvgSpec[band];

        nodeIn->m_Attributes.m_AvgSpec[band] = (area_1*avgSpec_1 + area_2*avgSpec_2)/(area_1 + area_2);
    }
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

    /**TODO Maybe not usefull to use sqrt. We do not need the rea lvalue, only the square value is necesary to sort cost*/
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

template< typename TGraph >
bool
SmallRegionsMergingFilter<TGraph>
::CheckMutuality(NodeType* node_1, NodeType* node_2)
{
	/**Go through the edge*/
	//Loop over edge
	auto outputGraph = const_cast< GraphType * > (this->GetOutput() );
	NodeType* bestNode = nullptr;
	float minDistance = std::numeric_limits<float>::max();

	for(auto edgeIt = node_2->Begin(); edgeIt != node_2->End(); edgeIt++)
	{
		auto adjNode = outputGraph->GetNodeAt(edgeIt->m_TargetId);
		if(adjNode->m_Attributes.m_Area <= m_MinimalSurface)
		{
			//std::cout << "Noeud adjacent a une surface " << adjNode->m_Attributes.m_Area << std::endl;
			//Compute distance
			float curDistance = ComputeSpectralDistance(node_2, adjNode);

			//Update minDistance
			if(curDistance < minDistance)
			{
				//std::cout << "Update distance" << std::endl;
				//Update min distance
				minDistance = curDistance;
				bestNode = adjNode;
			}else if(minDistance == curDistance)
			{
				//Get the node with the lowest ID
				if(adjNode->GetFirstPixelCoords() < bestNode->GetFirstPixelCoords())
				{
					//Update bestNode
					bestNode = adjNode;
				}
			}
		}
	}

	if(bestNode != nullptr)
	{
		if(bestNode->m_Id == node_1->m_Id)
		{
			return true;
		}
	}

	return false;

}
}//End obia
}//End otb

#endif
