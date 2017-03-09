#ifndef otbObiaBaatzSegmentationFilter_h
#define otbObiaBaatzSegmentationFilter_h
#include "itkImageRegionConstIterator.h"
#include <limits>

#include "otbObiaGraph.h"
#include "otbObiaGraphToGraphFilter.h"

namespace otb
{
namespace obia
{

/**Class specializing the merging cost function required by the generic filter*/
template< typename TCost, typename TGraph >
class BaatzMergingCost
{
public:
	/** Some convenient alias */
	using ValueType 	   = TCost;
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

	static ValueType Max(){return std::numeric_limits<TCost>::max();};

	/** Get/Set methods */
	void SetSpectralWeight(float spectralWeight){ m_SpectralWeight = spectralWeight;};
	void SetShapeWeight(float shapeWeight){ m_ShapeWeight = shapeWeight;};
	void SetBandWeights(std::vector<float> bandWeights){ m_BandWeights = bandWeights;};
	void SetThreshold(float threshold){m_Threshold = threshold;};

	const ValueType GetMax(){return m_MaxCost;};
	bool ComputeMergingCostsForThisNode(NodeType* curNode);
	bool ComputeMergingCostsForThisAdjNode(NodeType* curNode);
	ValueType ComputeMergingCost(NodeType* NodeIn, NodeType* NodeOut);
protected:

	// Maximal cost of merging
	ValueType m_MaxCost;

	// Relative inmportance given to the spectral information
	float m_SpectralWeight;

	// Relative importance given to geometrical information
	float m_ShapeWeight;

	// SThreshold
	float m_Threshold;

	// Relative importance weights for each band
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

	NodeType* GetBestAdjacentNode(NodeType* curNode);

	/** Get/Set methods */
	void SetGraph(GraphPointerType graph){ m_Graph = graph;};
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

	/**Update attributes of nodeIn*/
	void UpdateAttributes(NodeType * nodeIn, NodeType *  nodeOut);
};

template< typename TGraph >
class BaatzSegmentationFilter : public GraphToGraphFilter<TGraph, TGraph> 
{
public:

	/** Some convenient alias */
	using GraphType        = TGraph;
	using GraphPointerType = typename GraphType::Pointer;
	using NodeType         = typename GraphType::NodeType;
	using EdgeType         = typename GraphType::EdgeType;

	/** Standard class alias */
	using Self         = BaatzSegmentationFilter;
	using Superclass   = GraphToGraphFilter<TGraph, TGraph>;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self>;

	/** Method for creation through the object factory. */
  	itkNewMacro(Self);

  	/** Run-time type information (and related methods). */
  	itkTypeMacro(BaatzSegmentationFilter, GraphToGraphFilter);

  	void SetBandWeights(const std::vector<float>& weights){ m_BandWeights = weights; }
  	itkSetMacro(MaxNumberOfIterations, uint32_t);
  	itkGetMacro(SegmentationOver, bool);
  	itkSetMacro(Threshold, float);
  	itkSetMacro(SpectralWeight, float);
  	itkSetMacro(ShapeWeight, float);

protected:

	BaatzSegmentationFilter();

	virtual ~BaatzSegmentationFilter();

	void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

	void GenerateData();

private:

	void UpdateMergingCosts();
	float ComputeMergingCost(NodeType* n1, NodeType* n2);
	bool DoOneIteration();
	NodeType * GetBestAdjacentNode(NodeType * node);
	void Merge(NodeType* nodeIn, NodeType * nodeOut);
	void UpdateSpecificAttributes(NodeType * nodeIn,NodeType * nodeOut);

	// The maximum number of iterations to apply
	uint32_t m_MaxNumberOfIterations;

	// Flag indicating id the segmentation process finishes because
	// the segmentation is over or because we reach the maximum
	// number of iterations.
	bool m_SegmentationOver;

	// Threshold to determine if two adjacent nodes are similar
	float m_Threshold;

	// Relative inmportance given to the spectral information
	float m_SpectralWeight;

	// Relative importance given to geometrical information
	float m_ShapeWeight;

	// Relative importance weights for each band
	std::vector<float> m_BandWeights;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaBaatzSegmentationFilter.txx"
#endif
