#ifndef otbObiaLSSmallRegionsMergingFilter_h
#define otbObiaLSSmallRegionsMergingFilter_h
#include <limits>

#include "otbObiaLSBaatzSegmentationScheduler.h"
#include "otbObiaLSGraphToGraphFilter.h"
#include "itkMacro.h"
#include "itkSmartPointer.h"
//#include "otbObiaSmallRegionsMergingGraph.h"

namespace otb
{
namespace obia
{

/** \class SmallRegionMergingFilter
 *	\brief Class that builds an adjacency graph where all smalls regions are merged
 *
 */
//Define the graph type : SRM (Small Region Merging)
//using GraphType = Graph< Node< SRMNodeAttribute, SRMEdgeAttribute > >;

template< typename TGraph >
class LSSmallRegionsMergingScheduler : public LSGraphToGraphFilter< TGraph, TGraph >
{

public:
	/** Standard class alias */
	using Self 		   = LSSmallRegionsMergingScheduler;
	using SuperClass   = itk::Object;
	using Pointer      = itk::SmartPointer< Self >;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Some convenient alias */
	using GraphType                   = TGraph;
	using NodeType                    = typename GraphType::NodeType;
	using GraphPointerType            = typename GraphType::Pointer;
	using GraphOperationsType         = GraphOperations<GraphType>;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(LSSmallRegionsMergingScheduler, itk::LightObject);

	/** Get/Set methods */
	void SetMinimalSurface(const uint32_t minimalSurface){m_MinimalSurface = minimalSurface;};
	void SetNumberOfIterations(const uint32_t numberOfIterations){m_NumberOfIterations = numberOfIterations;};
protected:

	/** Constructor */
	LSSmallRegionsMergingScheduler();

	/** Destructor */
	~LSSmallRegionsMergingScheduler();

	/** Generation method */
	virtual void GenerateData();

	virtual void ComputePaddingValue();

	virtual void NoTilingExecution();

	virtual void TilingExecution();

private:

	enum FusionState{
		PARTIAL_FUSION = 0, // Next step is a partial fusion
		WAIT, // When all small regions of a graph merge, wait others graph
		AGGR_AND_FUSION, // Next step is merging all the graphs and achieve the fusion
		END // Next step is merging the graphs.
	};


	/** Extract stability Margins*/
	void ExtractStabilityMargins();

	/** Share Stability Margins*/
	std::vector< char > ShareStabilityMargins();

	/** Aggregate margins */
	void AggregateStabilityMargins();

	/**Merging small regions*/
	int PartialFusion(unsigned long int& accumulatedMemory);

	/** Extract small regions node*/
	std::vector< NodeType* > ExtractSmallRegionsNodes(std::vector< NodeType* > borderNodes);

	/**Compute max depth*/
	uint32_t ComputeMaxDepth();

	/**Recurisve Exploration Depth*/
	void RecursiveSmallRegionsDepthBreadthFirstExploration(NodeType* node,
														   std::unordered_map< NodeType*, uint32_t >& visitedNodes,
														   const GraphPointerType graph,
														   const uint32_t currentNumberOfAdjacencyLayers,
														   const uint32_t currentArea);

	/** Check if node is in border*/
	bool IsBorderNode(NodeType* node, const std::vector< typename GraphOperations<TGraph>::NodeType* > borderNodes);

	/** Check if duplicated node*/
	bool HasDuplicatedNodes();

	/**Create output*/
	void CreateOutput();

	/** Convert Graph to Image*/
	void ConvertToImage(uint32_t ty, uint32_t tx);

	/**Surface value to consider a small region*/
	uint32_t m_MinimalSurface;

	/** Internal member attributes */

	// The serialized stability margin extracted from the current graph
	std::vector< char > m_SerializedStabilityMargin;

	// The maximum number of bytes to send via MPI
	unsigned long int m_MaxNumberOfBytes;

	// The number of iteration
	int m_NumberOfIterations;

	// Flag indicating if fusion is over
	bool m_FusionOver;

};
}//End obia
}//End otb
#include "otbObiaLSSmallRegionsMergingScheduler.txx"
#endif
