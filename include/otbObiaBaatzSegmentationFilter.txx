#ifndef otbObiaBaatzSegmentationFilter_txx
#define otbObiaBaatzSegmentationFilter_txx
#include "otbObiaBaatzSegmentationFilter.h"

namespace otb
{
namespace obia
{

template< typename TCost, typename TGraph >
BaatzMergingCost<TCost, TGraph>
::BaatzMergingCost()
{

}

template< typename TCost, typename TGraph >
BaatzMergingCost<TCost, TGraph>
::~BaatzMergingCost()
{
}

template< typename TCost, typename TGraph >
bool
BaatzMergingCost<TCost, TGraph>::ComputeMergingCostsForThisNode(NodeType* curNode)
{
    (void) curNode;
    return true;
}

template< typename TCost, typename TGraph >
bool
BaatzMergingCost<TCost, TGraph>::ComputeMergingCostsForThisAdjNode(NodeType* curNode)
{
    (void) curNode;
    return true;
}

template< typename TCost, typename TGraph >
typename BaatzMergingCost< TCost, TGraph>::ValueType
BaatzMergingCost<TCost, TGraph>::ComputeMergingCost(NodeType* n1, NodeType* n2)
{
    // Retrieve attributes
    const float aArea = n1->m_Attributes.m_Area;
    const float bArea = n2->m_Attributes.m_Area;
    const float areaSum = aArea + bArea;


    float mean;
    float stddev;
    float stddevTmp;
    float colorF;
    float colorH = 0.0f;

    if(m_BandWeights.size() < 1)
    {
        m_BandWeights.assign(n1->m_Attributes.m_AvgSpec.size(), 1.0f);
    }

    for (uint32_t band = 0; band < n1->m_Attributes.m_AvgSpec.size(); band++)
    {
        mean = (n1->m_Attributes.m_AvgSpec[band] * aArea + n2->m_Attributes.m_AvgSpec[band] * bArea) / areaSum;
        stddevTmp = aArea * n1->m_Attributes.m_StdSpec[band] * n1->m_Attributes.m_StdSpec[band];
        stddevTmp += bArea * n2->m_Attributes.m_StdSpec[band] * n2->m_Attributes.m_StdSpec[band];
        stddevTmp += aArea * (n1->m_Attributes.m_AvgSpec[band] - mean) * (n1->m_Attributes.m_AvgSpec[band] - mean);
        stddevTmp += bArea * (n2->m_Attributes.m_AvgSpec[band] - mean) * (n2->m_Attributes.m_AvgSpec[band] - mean);

        stddev = std::sqrt(stddevTmp / areaSum);
        colorF = n1->m_Attributes.m_Area * n1->m_Attributes.m_StdSpec[band] + n2->m_Attributes.m_Area * n2->m_Attributes.m_StdSpec[band];

        colorF = m_BandWeights[band] * ((areaSum * stddev) - colorF);
        colorH += colorF;

    }

    float H = m_SpectralWeight * colorH;

    if(H < this->m_Threshold)
    {
        const float aPerimeter = n1->m_Attributes.m_Perimeter;
        const float bPerimeter = n2->m_Attributes.m_Perimeter;

        // Compute perimeter of the merged nodes
        const float boundaryLength = (n1->FindEdge(n2->m_Id))->m_Boundary;
        const float mergedPerimeter = aPerimeter + bPerimeter - 2 * boundaryLength;

        // Compute bounding box lengths
        const float aBboxLen = 2*(n1->m_BoundingBox[2] + n1->m_BoundingBox[3]);
        const float bBboxLen = 2*(n2->m_BoundingBox[2] + n2->m_BoundingBox[3]);

        // Merged bounding box
        auto mergedBbox = SpatialTools::GetMergedBoundingBox(n1->m_BoundingBox, n2->m_BoundingBox);
        const float mergedBboxLen = 2*(mergedBbox[2] + mergedBbox[3]);

        // Smoothness factor
        const float smoothF = ( (areaSum * mergedPerimeter) / mergedBboxLen) - ( ( (aArea * aPerimeter) / aBboxLen) + ( (bArea * bPerimeter) / bBboxLen));

        // Compactness factor
        const float compactF = areaSum * mergedPerimeter / std::sqrt(areaSum) - ( aArea * aPerimeter / std::sqrt(aArea) +
                               bArea * bPerimeter / std::sqrt(bArea) );

        // Spatial heterogeneity
        float spatialH = m_ShapeWeight * compactF + (1 - m_ShapeWeight) * smoothF;

        return (H + (1-this->m_SpectralWeight)*spatialH);
    }
    else
    {
        return H;
    }
}

template< typename TGraph >
BaatzHeuristic<TGraph>
::BaatzHeuristic()
{

}

template< typename TGraph >
BaatzHeuristic<TGraph>
::~BaatzHeuristic()
{

}
template< typename TGraph >
typename BaatzHeuristic<TGraph>::NodeType*
BaatzHeuristic<TGraph>
::GetBestAdjacentNode(NodeType* node)
{

    auto outputGraph = this->m_Graph;

    if(node->m_Valid)
    {
        auto cost = node->m_Edges.front().m_Attributes.m_MergingCost;

        if(cost < this->m_Threshold)
        {
            auto bestAdjNode = outputGraph->GetNodeAt(node->m_Edges.front().m_TargetId);

            if(bestAdjNode->m_Valid)
            {
                auto bbNode = outputGraph->GetNodeAt(bestAdjNode->m_Edges.front().m_TargetId);

                if(bbNode->m_Id == node->m_Id)
                {
                    if(bestAdjNode->GetFirstPixelCoords() < node->GetFirstPixelCoords())
                    {
                        return bestAdjNode;
                    }
                    else
                    {
                        return node;
                    }
                }
                else
                {
                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }
        else
        {
            return nullptr;
        }
    }
    else
    {
        return nullptr;
    }
}


template< typename TGraph >
BaatzUpdateAttribute<TGraph>
::BaatzUpdateAttribute()
{
}

template< typename TGraph >
BaatzUpdateAttribute<TGraph>
::~BaatzUpdateAttribute()
{
}

template<typename TGraph >
void
BaatzUpdateAttribute<TGraph>
::UpdateAttributes(NodeType * nodeIn, NodeType *  nodeOut)
{
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
BaatzSegmentationFilter<TGraph>
::BaatzSegmentationFilter() : 
m_MaxNumberOfIterations(75), m_SegmentationOver(false), m_Threshold(1600), m_SpectralWeight(0.5), m_ShapeWeight(0.5)
{
    m_BandWeights.clear();
}

template< typename TGraph >
BaatzSegmentationFilter<TGraph>
::~BaatzSegmentationFilter()
{
}

template< typename TGraph >
void 
BaatzSegmentationFilter<TGraph>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

template< typename TGraph >
void 
BaatzSegmentationFilter<TGraph>
::GenerateData()
{
    // Graft the input to the output
    // #TODO it may be good to create an abstract filter otbObiaInPlaceGraphFilter for this purpose
    // by following the same scheme of itkInPlaceImageFilter.
    GraphType * inputGraphPtr = dynamic_cast< GraphType * >( this->itk::ProcessObject::GetInput(0));

    auto outputGraph = this->GetOutput();
    outputGraph->GraftGraphByMove(inputGraphPtr);
    // Release the input graph.
    inputGraphPtr->Reset();

    // Reconditionning of the graph
    for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {
        nodeIt->m_HasToBeRemoved = false;
        nodeIt->m_Valid = true;
        nodeIt->m_Attributes.m_HasPreviouslyMerged = true;
        for(auto& edg : nodeIt->m_Edges)
        {
            edg.m_Attributes.m_MergingCost = std::numeric_limits<float>::max();
            edg.m_Attributes.m_CostUpdated = false;
        }
    }

    if(m_BandWeights.size() < 1)
    {
        m_BandWeights.assign(inputGraphPtr->GetNumberOfSpectralBands(), 1.0f);
    }

    for(uint32_t i = 0; i < m_MaxNumberOfIterations; i++)
    {
        std::cout << "iteration " << i+1 << std::endl;
        if(!DoOneIteration())
        {
            m_SegmentationOver = true;
            break;
        }
    }
}

template< typename TGraph >
bool
BaatzSegmentationFilter<TGraph>
::DoOneIteration()
{
    auto outputGraph = this->GetOutput();
    bool merged = false;

    /* Update the costs of merging between adjacent nodes */
    UpdateMergingCosts();

    for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {
        auto nodeIn = GetBestAdjacentNode(&(*nodeIt));

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

template< typename TGraph >
void
BaatzSegmentationFilter<TGraph>
::Merge(NodeType* nodeIn, NodeType * nodeOut)
{
    UpdateSpecificAttributes(nodeIn, nodeOut);
    auto outputGraph = this->GetOutput();
    outputGraph->Merge(nodeIn, nodeOut);
}

template< typename TGraph >
void
BaatzSegmentationFilter<TGraph>
::UpdateSpecificAttributes(NodeType * nodeIn, NodeType * nodeOut)
{
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

template<typename TGraph >
typename BaatzSegmentationFilter<TGraph>::NodeType*
BaatzSegmentationFilter<TGraph>
::GetBestAdjacentNode(NodeType * node)
{
    auto outputGraph = this->GetOutput();

    if(node->m_Valid)
    {
        auto cost = node->m_Edges.front().m_Attributes.m_MergingCost;

        if(cost < this->m_Threshold)
        {
            auto bestAdjNode = outputGraph->GetNodeAt(node->m_Edges.front().m_TargetId);

            if(bestAdjNode->m_Valid)
            {
                auto bbNode = outputGraph->GetNodeAt(bestAdjNode->m_Edges.front().m_TargetId);

                if(bbNode->m_Id == node->m_Id)
                {
                    if(bestAdjNode->GetFirstPixelCoords() < node->GetFirstPixelCoords())
                    {
                        return bestAdjNode;
                    }
                    else
                    {
                        return node;
                    }
                }
                else
                {
                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }
        else
        {
            return nullptr;
        }
    }
    else
    {
        return nullptr;
    }
}

template< typename TGraph >
void
BaatzSegmentationFilter<TGraph>
::UpdateMergingCosts()
{
    auto outputGraph = this->GetOutput();

    float min_cost;
    uint64_t min_id  = 0;
    std::size_t idx, min_idx;

    for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {
        for(auto edgeIt = nodeIt->m_Edges.begin(); edgeIt != nodeIt->m_Edges.end(); edgeIt++)
        {
            edgeIt->m_Attributes.m_CostUpdated = false;
        }
    }

    for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {

        min_cost = std::numeric_limits<float>::max();
        idx = 0;
        min_idx = 0;
        nodeIt->m_HasToBeRemoved = false;
        nodeIt->m_Valid = true;

        for(auto edgeIt = nodeIt->m_Edges.begin(); edgeIt != nodeIt->m_Edges.end(); edgeIt++)
        {
            auto adjNode = outputGraph->GetNodeAt(edgeIt->m_TargetId);

            if(!(edgeIt->m_Attributes.m_CostUpdated) && 
              (nodeIt->m_Attributes.m_HasPreviouslyMerged || 
               adjNode->m_Attributes.m_HasPreviouslyMerged))
            {
                // Compute the cost.
                edgeIt->m_Attributes.m_MergingCost = ComputeMergingCost(&(*nodeIt), adjNode);
                auto adjNodeToCurr = adjNode->FindEdge(nodeIt->m_Id);
                adjNodeToCurr->m_Attributes.m_MergingCost = edgeIt->m_Attributes.m_MergingCost;
                edgeIt->m_Attributes.m_CostUpdated = true;
                adjNodeToCurr->m_Attributes.m_CostUpdated = true;
            }


            if(edgeIt->m_Attributes.m_MergingCost < min_cost)
            {
                min_cost = edgeIt->m_Attributes.m_MergingCost;
                min_id = adjNode->GetFirstPixelCoords();
                min_idx = idx;
            }
            else if(min_cost == edgeIt->m_Attributes.m_MergingCost)
            {
                if(adjNode->GetFirstPixelCoords() < min_id)
                {
                    min_id = adjNode->GetFirstPixelCoords();
                    min_idx = idx;
                }
            }

            idx++;    
        }

        std::swap(nodeIt->m_Edges[0], nodeIt->m_Edges[min_idx]);
                
    }

    // Reset the merge flag for all the regions.
    for(auto nodeIt =outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    {
        nodeIt->m_Attributes.m_HasPreviouslyMerged = false;
    }
}

template< typename TGraph >
float
BaatzSegmentationFilter<TGraph>
::ComputeMergingCost(NodeType* n1, NodeType* n2)
{
    // Retrieve attributes
    const float aArea = n1->m_Attributes.m_Area;
    const float bArea = n2->m_Attributes.m_Area;
    const float areaSum = aArea + bArea;


    float mean;
    float stddev;
    float stddevTmp;
    float colorF;
    float colorH = 0.0f;

    for (uint32_t band = 0; band < n1->m_Attributes.m_AvgSpec.size(); band++)
    {
        mean = (n1->m_Attributes.m_AvgSpec[band] * aArea + n2->m_Attributes.m_AvgSpec[band] * bArea) / areaSum;
        stddevTmp = aArea * n1->m_Attributes.m_StdSpec[band] * n1->m_Attributes.m_StdSpec[band];
        stddevTmp += bArea * n2->m_Attributes.m_StdSpec[band] * n2->m_Attributes.m_StdSpec[band];
        stddevTmp += aArea * (n1->m_Attributes.m_AvgSpec[band] - mean) * (n1->m_Attributes.m_AvgSpec[band] - mean);
        stddevTmp += bArea * (n2->m_Attributes.m_AvgSpec[band] - mean) * (n2->m_Attributes.m_AvgSpec[band] - mean);
        stddev = std::sqrt(stddevTmp / areaSum);
        colorF = n1->m_Attributes.m_Area * n1->m_Attributes.m_StdSpec[band] + n2->m_Attributes.m_Area * n2->m_Attributes.m_StdSpec[band];
        colorF = m_BandWeights[band] * ((areaSum * stddev) - colorF);
        colorH += colorF;
    }

    float H = m_SpectralWeight * colorH;


    if(H < this->m_Threshold)
    {
        const float aPerimeter = n1->m_Attributes.m_Perimeter;
        const float bPerimeter = n2->m_Attributes.m_Perimeter;

        // Compute perimeter of the merged nodes
        const float boundaryLength = (n1->FindEdge(n2->m_Id))->m_Boundary;
        const float mergedPerimeter = aPerimeter + bPerimeter - 2 * boundaryLength;

        // Compute bounding box lengths
        const float aBboxLen = 2*(n1->m_BoundingBox[2] + n1->m_BoundingBox[3]);
        const float bBboxLen = 2*(n2->m_BoundingBox[2] + n2->m_BoundingBox[3]);

        // Merged bounding box
        auto mergedBbox = SpatialTools::GetMergedBoundingBox(n1->m_BoundingBox, n2->m_BoundingBox);
        const float mergedBboxLen = 2*(mergedBbox[2] + mergedBbox[3]);
      
        // Smoothness factor
        const float smoothF = ( (areaSum * mergedPerimeter) / mergedBboxLen) - ( ( (aArea * aPerimeter) / aBboxLen) + ( (bArea * bPerimeter) / bBboxLen));

        // Compactness factor
        const float compactF = areaSum * mergedPerimeter / std::sqrt(areaSum) - ( aArea * aPerimeter / std::sqrt(aArea) + 
                               bArea * bPerimeter / std::sqrt(bArea) );
      
        // Spatial heterogeneity
        float spatialH = m_ShapeWeight * compactF + (1 - m_ShapeWeight) * smoothF;

        return (H + (1-this->m_SpectralWeight)*spatialH);
    }
    else
    {
        return H;
    }
}

}
}


#endif
