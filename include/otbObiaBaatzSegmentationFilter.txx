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
BaatzMergingCost<TCost, TGraph>::ComputeMergingCostsForThisNode(const NodeType* curNode)
{
    return true;
}

template< typename TCost, typename TGraph >
bool
BaatzMergingCost<TCost, TGraph>::ComputeMergingCostsForThisAdjNode(const NodeType* curNode)
{
    return true;
}

template< typename TCost, typename TGraph >
typename BaatzMergingCost< TCost, TGraph>::ValueType
BaatzMergingCost<TCost, TGraph>::ComputeMergingCost(NodeType* n1, NodeType* n2)  // TODO: should be const
{
    // Retrieve attributes
    const float aArea = n1->m_Attributes.m_Area;
    const float bArea = n2->m_Attributes.m_Area;
    const float areaSum = aArea + bArea;
    
    float colorH = 0.0f;

    // TODO: No. We should do this elsewhere!
   //    if(m_BandWeights.size() < 1)
   //    {
   //        m_BandWeights.assign(n1->m_Attributes.m_AvgSpec.size(), 1.0f);
   //    }

    for (uint32_t band = 0; band < n1->m_Attributes.m_AvgSpec.size(); band++)
    {
        float mean;
        float stddev;
        float stddevTmp;
        float colorF;

        mean = (n1->m_Attributes.m_AvgSpec[band] * aArea + n2->m_Attributes.m_AvgSpec[band] * bArea) / areaSum;
        stddevTmp = aArea * n1->m_Attributes.m_StdSpec[band] * n1->m_Attributes.m_StdSpec[band];
        stddevTmp += bArea * n2->m_Attributes.m_StdSpec[band] * n2->m_Attributes.m_StdSpec[band];
        stddevTmp += aArea * (n1->m_Attributes.m_AvgSpec[band] - mean) * (n1->m_Attributes.m_AvgSpec[band] - mean);
        stddevTmp += bArea * (n2->m_Attributes.m_AvgSpec[band] - mean) * (n2->m_Attributes.m_AvgSpec[band] - mean);

        stddev = std::sqrt(stddevTmp / areaSum);
        colorF = n1->m_Attributes.m_Area * n1->m_Attributes.m_StdSpec[band] + n2->m_Attributes.m_Area * n2->m_Attributes.m_StdSpec[band];

        // TODO: replace it, when m_BandWeights will be properly set elsewhere
//     colorF = m_BandWeights[band] * ((areaSum * stddev) - colorF);
       colorF = ((areaSum * stddev) - colorF);




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

  if(node->m_Valid
      && node->m_Edges.size() > 0) // Since the introducing of no-data, a node can be alone
  {
    auto cost = node->m_Edges.front().m_Attributes.m_MergingCost;

    if(cost < this->m_Threshold)
    {
      auto bestAdjNode = outputGraph->GetNodeAt(node->m_Edges.front().m_TargetId);

      if(bestAdjNode->m_Valid
          && bestAdjNode->m_Edges.size() > 0) // Since the introducing of no-data, a node can be alone
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
      }
    }
  }
  return nullptr;

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

    for(unsigned int b = 0; b < nodeIn->m_Attributes.m_AvgSpec.size(); ++b)
    {
        float mean;
        float stddev;
        float stddevTmp;

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
}
}

#endif
