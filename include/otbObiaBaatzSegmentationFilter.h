#ifndef otbObiaBaatzSegmentationFilter_h
#define otbObiaBaatzSegmentationFilter_h
#include "itkImageRegionConstIterator.h"
#include <limits>

#include "otbObiaGraph.h"
#include "otbObiaGraphToGraphFilter.h"

/**
\file otbObiaBaatzSegmentationFilter.h
\brief This file define all the classes (merging cost, heuristic, update)
	   required in order to compute Baatz & Schäpe Segmentation.
*/
namespace otb
{
namespace obia
{

/**\class BaatzMergingCost otbObiaBaatzSegmentationFilter.h
 * \brief Class specializing the merging cost function required by the generic filter*/
template< typename TCost, typename TGraph >
class BaatzMergingCost
{
public:
    /** Some convenient alias */
    using ValueType        = TCost;
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = BaatzMergingCost;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    BaatzMergingCost();
    ~BaatzMergingCost();

    /**\brief Return the max value according to the input template type: TCost.
     * For example, it can return the max double value if TCost = double*/
    static ValueType Max(){return std::numeric_limits<TCost>::max();};

    /**\brief Set the attribut m_SpectralWeight */
    void SetSpectralWeight(float spectralWeight){ m_SpectralWeight = spectralWeight;};

    /**\brief Set the attribut m_ShapeWeight */
    void SetShapeWeight(float shapeWeight){ m_ShapeWeight = shapeWeight;};

    /**\brief Set the attribut m_BandWeights */
    void SetBandWeights(std::vector<float> bandWeights){ m_BandWeights = bandWeights;};

    /**\brief Set the attribut m_Threshold */
    void SetThreshold(float threshold){m_Threshold = threshold;};

    /**\brief Return the max cost */
    const ValueType GetMax() const {return m_MaxCost;};

    /**\brief This method indicate if it is necessary to compute merging cost for the given node
     * @param: Node to test
     * @return : return true if the merging cost is required*/
    bool ComputeMergingCostsForThisNode(const NodeType* curNode);

    /**\brief This method indicate if it is necessary to compute merging cost for the adjacent node
     * @param: Node to test
     * @return : return true if the merging cost is required*/
    bool ComputeMergingCostsForThisAdjNode(const NodeType* curNode);

    /**\fn ValueType ComputeMergingCost(NodeType* NodeIn, NodeType* NodeOut);
     * \brief This method compute the merging cost for Baatz & Schäpe segmentation.
     * The cost is computed in this way:
     * - Compute radiometric criterion:
     * 		-# For each band \n
     * 		  	-# compute mean of final node (mean1*area1 + mean2*area2)/(area1 + area2)
     * 			-# compute standard deviation
     * 			-# compute area1*stdev1 + area2*stdev2
     * 			-# compute weighted value : colorH += m_BandWeights[band] * ((areaSum * stddev) - colorF) \n
     * 		-# Final radiometric criterion = m_SpectralWeight*colorH \n
     * - Compute geometric criterion:
     * 		-# For each band \n
     * 			-# Compute smoothness factor smoothF
     * 			-# Compute compactness factor compactF
     * 			-# Compute weighted value : m_ShapeWeight[band]*compacF + (1 - m_ShapeWeight[band])*smoothF \n
     * \param: NodeIn
     * \param: NodeOut
     * \return : return true if the merging cost is required*/
    ValueType ComputeMergingCost(NodeType* NodeIn, NodeType* NodeOut); // TODO: should be const
protected:

    /**\brief Maximal cost of merging*/
    ValueType m_MaxCost;

    /**\brief Relative inmportance given to the spectral information*/
    float m_SpectralWeight;

    /**\Relative importance given to geometrical information*/
    float m_ShapeWeight;

    /**\ Threshold for computing merging cost */
    float m_Threshold;

    /**\brief Relative importance weights for each band */
    std::vector<float> m_BandWeights;

};

/**Class specializing the heuristic function required by the generic filter*/
template<typename TGraph >
class BaatzHeuristic
{

public:
    /** Some convenient alias */
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = BaatzHeuristic;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;


    /** Method for creation through the object factory. */
    //itkNewMacro(Self);
    BaatzHeuristic();
    ~BaatzHeuristic();

    /**\brief Return the best adjacent node (minimal cost)
     * @param : Node
     * @return : best adjacent node*/
    NodeType* GetBestAdjacentNode(NodeType* curNode);

    /**\brief Set the graph pointer
     * \param : graph pointer*/
    void SetGraph(GraphPointerType graph){ m_Graph = graph;};

    /**\brief Set the threshold
     * \param : threshold*/
    void SetThreshold(float threshold){m_Threshold = threshold;};

protected:

    /** Pointer to the graph*/
    GraphPointerType m_Graph;

    /** Threshold for Baatz decision*/
    float m_Threshold;
};

template<typename TGraph >
class BaatzUpdateAttribute
{

public:
    /** Some convenient alias */
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = BaatzUpdateAttribute;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    BaatzUpdateAttribute();
    ~BaatzUpdateAttribute();

    /**\brief Update attributes of nodeIn using NodeOut information
     * For example, it will update the area value (areaSum = areaNodeIn + areaNodeOut) and also
     * the mean value (using (meanValue1*area1 + meanValue2*area2)/areaSum)\n
     * @param: Node In
     * @param: Node Out (which merged with nodeIn)*/
    void UpdateAttributes(NodeType * nodeIn, NodeType *  nodeOut);
};


} // end of namespace obia
} // end of namespace otb
#include "otbObiaBaatzSegmentationFilter.txx"
#endif
