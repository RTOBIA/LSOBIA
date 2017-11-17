#ifndef __otbObiaGenericRegionMergingFilter_txx
#define __otbObiaGenericRegionMergingFilter_txx
#include "otbObiaGenericRegionMergingFilter.h"
#include "itkTimeProbe.h"
#include <thread>
#include <iostream>
#include <vector>
namespace otb
{
namespace obia
{

/*
 * Class to help with threads.
 * We set the range of elements to process within the array.
 * This can process:
 * -nodes
 * -pairs of nodes
 */
template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
class ThreadWorker
{
public:

  ThreadWorker(){};
  ~ThreadWorker(){};

  void SetGraph(TOutputGraph* graph)
  {
    m_OutputGraph = graph;
  }
  void SetMergingCostFunc(TMergingCostFunc * mergingCostFunc)
  {
    m_MergingCostFunc = mergingCostFunc;
  }
  void SetUpdateAttributeFunc(TUpdateAttributeFunc * func)
  {
    m_UpdateAttributesFunc = func;
  }
  void SetRange(uint64_t istart, uint64_t iend)
  {
    m_RangeStart = istart;
    m_RangeEnd = iend;
  }
  void SetPairsOfNodesToProcess(std::vector<std::pair<typename TOutputGraph::NodeType * ,typename TOutputGraph::NodeType * >> &pairs)
  {
    m_PairsOfNodesToProcess = pairs;
  }

  /*
   * Compute merging costs of the nodes
   */
  void ComputeMergingCosts()
  {
    uint64_t minNodeId = m_OutputGraph->GetNumberOfNodes()+1, idx, minIdx;
    typename TMergingCostFunc::ValueType minCost;

    // Loop over the nodes
    for(auto nodeIt = m_OutputGraph->Begin() + m_RangeStart; nodeIt != m_OutputGraph->Begin() + m_RangeEnd; nodeIt++)
      {
      nodeIt->m_ThreadSafe = true;
      nodeIt->m_ThreadSafeForMerge = true;

      if(m_MergingCostFunc->ComputeMergingCostsForThisNode(&(*nodeIt))
          && nodeIt->m_Edges.size()>0) // Since the introducing of no-data, a node can be alone
        {

        // The merging cost function must give the maximum value of the merging cost.
        minCost = TMergingCostFunc::Max();
        idx = 0;
        minIdx = 0;
        nodeIt->m_HasToBeRemoved = false;
        nodeIt->m_Valid = true;

        // Loop over the edges
        for(auto& edgeIt : nodeIt->m_Edges)
          {

          // Tell if this node has edges targeting nodes which are outside the thread range
          if ( edgeIt.m_TargetId < m_RangeStart || edgeIt.m_TargetId >= m_RangeEnd )
            {
            nodeIt->m_ThreadSafe = false;
            nodeIt->m_ThreadSafeForMerge = false;
            }

          // Retrieve the adjacent node.
          auto adjNode = m_OutputGraph->GetNodeAt(edgeIt.m_TargetId);
          if(m_MergingCostFunc->ComputeMergingCostsForThisAdjNode(adjNode))
            {

            // If the cost is not updated and if one of the adjacent nodes
            // has merged at the previous iteration then we must compute the
            // merging cost.
            if(nodeIt->m_Attributes.m_HasPreviouslyMerged || adjNode->m_Attributes.m_HasPreviouslyMerged)
              {
              edgeIt.m_Attributes.m_MergingCost = m_MergingCostFunc->ComputeMergingCost(&(*nodeIt), adjNode);
              }

            // If the current cost is minimum than we record it.
            if(edgeIt.m_Attributes.m_MergingCost < minCost)
              {
              minCost = edgeIt.m_Attributes.m_MergingCost;
              minNodeId = adjNode->GetFirstPixelCoords();
              minIdx = idx;
              }

            // In case of equality, we keep the adjacent node with the lower starting
            // coordinates.
            else if(minCost == edgeIt.m_Attributes.m_MergingCost)
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

    // grow the thread-unsafe-for-merge region (we know that the merge operation requires to R/W adjacent nodes of
    // the pair of nodes to merge)
    for(auto nodeIt = m_OutputGraph->Begin() + m_RangeStart; nodeIt != m_OutputGraph->Begin() + m_RangeEnd; nodeIt++)
      {
      if (nodeIt->m_ThreadSafe == false)
        {
        for(auto& edgeIt : nodeIt->m_Edges)
          {
          if (m_RangeStart <= edgeIt.m_TargetId && edgeIt.m_TargetId < m_RangeEnd)
            {
            auto adjNode = m_OutputGraph->GetNodeAt(edgeIt.m_TargetId);
            adjNode->m_ThreadSafeForMerge = false;
            }
          }
        }
      }


  } // ComputeMergingCosts()

  /*
   * Reset the merge flag(s) of nodes
   */
  void ResetMergeFlags()
  {
    for(auto nodeIt = m_OutputGraph->Begin() + m_RangeStart; nodeIt != m_OutputGraph->Begin() + m_RangeEnd; nodeIt++)
      {
      nodeIt->m_Attributes.m_HasPreviouslyMerged = false;
      }
  } // ResetMergeFlags()

  /*
   * Performs the merge of a pair of nodes
   */
  void FusionOfPairs()
  {
    for (auto pairIt =  m_PairsOfNodesToProcess.begin() + m_RangeStart ; pairIt != m_PairsOfNodesToProcess.begin() + m_RangeEnd; pairIt++)
      {
      typename TOutputGraph::NodeType * nodeIn = pairIt->first;
      typename TOutputGraph::NodeType * nodeOut = pairIt->second;

      // Update attributes of thread-safe nodes
      m_UpdateAttributesFunc->UpdateAttributes(nodeIn, nodeOut);

      // Fusion of the bounding box
      SpatialTools::MergeBoundingBox(nodeIn->m_BoundingBox, nodeOut->m_BoundingBox);

      // Merge the edges
      m_OutputGraph->MergeEdge(nodeIn, nodeOut);

      // Fusion of the contour
      nodeIn->m_Contour.MergeWith(nodeOut->m_Contour, m_OutputGraph->GetImageWidth(), m_OutputGraph->GetImageHeight());

      // nodeOut has to be removed, so we mark it as it
      nodeOut->m_HasToBeRemoved = true;
      }
  } // FusionOfPairs()

  /*
   * Reconditionning of the graph
   */
  void Reconditioning()
  {
    for(auto nodeIt = m_OutputGraph->Begin() + m_RangeStart; nodeIt != m_OutputGraph->Begin() + m_RangeEnd; nodeIt++)
    {
      nodeIt->m_HasToBeRemoved = false;
      nodeIt->m_Valid = true;
      nodeIt->m_Attributes.m_HasPreviouslyMerged = true;
      for(auto& edg : nodeIt->m_Edges)
      {
        edg.m_Attributes.m_MergingCost = m_MergingCostFunc->GetMax();
      }
    }
  }

private:
  TOutputGraph*           m_OutputGraph;
  TMergingCostFunc *      m_MergingCostFunc;
  TUpdateAttributeFunc *  m_UpdateAttributesFunc;
  int64_t                 m_RangeStart;
  int64_t                 m_RangeEnd;
  std::vector<std::pair<typename TOutputGraph::NodeType * ,typename TOutputGraph::NodeType * >>
                          m_PairsOfNodesToProcess;

};

template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
void ComputeMergingCostsInRange (TOutputGraph * graph, TMergingCostFunc * mergingCostFunc, uint64_t start, uint64_t end)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
  worker.SetMergingCostFunc( mergingCostFunc );
  worker.SetRange(start, end);
  worker.ComputeMergingCosts();

}

template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
void ResetMergeFlagsInRange (TOutputGraph * graph, TMergingCostFunc * mergingCostFunc, uint64_t start, uint64_t end)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
  worker.SetMergingCostFunc( mergingCostFunc );
  worker.SetRange(start, end);
  worker.ResetMergeFlags();

}

template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
void FusionOfPairsInRange (TOutputGraph * graph, TUpdateAttributeFunc * updateAttributeFunc,
    std::vector<std::pair<typename TOutputGraph::NodeType * ,typename TOutputGraph::NodeType * >> &pairsOfNodesToProcess,
    uint64_t start, uint64_t end)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
  worker.SetUpdateAttributeFunc( updateAttributeFunc );
  worker.SetPairsOfNodesToProcess(pairsOfNodesToProcess);
  worker.SetRange(start, end);
  worker.FusionOfPairs();

}

template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
void ReconditioningInRange (TOutputGraph * graph, TMergingCostFunc * mergingCostFunc, uint64_t start, uint64_t end)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
  worker.SetMergingCostFunc( mergingCostFunc );
  worker.SetRange(start, end);
  worker.Reconditioning();

}

template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc>
bool
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::DoOneIteration()
{
  // Output graph
  auto outputGraph = this->GetOutput();

////////////////////////////////////////////////////////////////////////////////
////            Check graph consistency.   TODO: remove this                 //
//  for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
//    {
//      assert(0 <= nodeIt->m_Id && nodeIt->m_Id < outputGraph->GetNumberOfNodes());
//      for (auto& edgeIt : nodeIt->m_Edges)
//        {
//          assert(0 <= edgeIt.m_TargetId && edgeIt.m_TargetId < outputGraph->GetNumberOfNodes());
//          auto adjacentNode = outputGraph->GetNodeAt(edgeIt.m_TargetId);
//          auto edgeFromAdjacentNodeToCurrentNode = adjacentNode->FindEdge(nodeIt->m_Id);
//          assert(edgeFromAdjacentNodeToCurrentNode != adjacentNode->m_Edges.end());
//        }
//    }
////////////////////////////////////////////////////////////////////////////////

  // Measure iteration total time
  itk::TimeProbe titeration;  titeration.Start();

  // Compute the merging costs for all the pairs of adjacent nodes.
  itk::TimeProbe tcosts;  tcosts.Start();
  ComputeMergingCosts();

//  //////////////////////////////////////////////////////////////////////////////
//  //            Check graph consistency.   TODO: remove this                 //
//  for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
//  {
//	  if (nodeIt->m_ThreadSafeForMerge)
//	  {
//		  for (auto& edgeIt : nodeIt->m_Edges)
//		  {
//			  auto adjacentNodeToThreadSafeForMergeNode = outputGraph->GetNodeAt(edgeIt.m_TargetId);
//			  assert(adjacentNodeToThreadSafeForMergeNode->m_ThreadSafe);
//		  }
//	  }
//	  assert(0 <= nodeIt->m_Id && nodeIt->m_Id < outputGraph->GetNumberOfNodes());
//	  for (auto& edgeIt : nodeIt->m_Edges)
//	  {
//		  assert(0 <= edgeIt.m_TargetId && edgeIt.m_TargetId < outputGraph->GetNumberOfNodes());
//		  auto adjacentNode = outputGraph->GetNodeAt(edgeIt.m_TargetId);
//		  auto edgeFromAdjacentNodeToCurrentNode = adjacentNode->FindEdge(nodeIt->m_Id);
//		  assert(edgeFromAdjacentNodeToCurrentNode != adjacentNode->m_Edges.end());
//	  }
//  }
//  //////////////////////////////////////////////////////////////////////////////

  tcosts.Stop(); timingsValues[1] += tcosts.GetTotal();
  std::cout << std::setprecision(3) << tcosts.GetTotal() << "\t";

  // Flag indicating if at least one merge has been done during the iteration.
  bool merged = false;

  itk::TimeProbe tmerge;  tmerge.Start();
  itk::TimeProbe tmergePre;  tmergePre.Start();

  // First part is not multithreaded, we just identify pairs of nodes to merge
  typedef typename TInputGraph::NodeType NodeType;
  typedef std::pair<NodeType * ,NodeType * > PairOfNodes;
  std::vector<PairOfNodes> pairsOfNodesToProcess_parallel;
  std::vector<PairOfNodes> pairsOfNodesToProcess_sequential;
  std::vector<PairOfNodes> pairsOfNodesToProcess_all;

  // Let's allocate the maximum number of pairs possible (ie 1+n/2)
  const int64_t maxNbOfPairs = outputGraph->GetNumberOfNodes() / 2 + 1;
  pairsOfNodesToProcess_parallel.reserve(maxNbOfPairs);
  pairsOfNodesToProcess_sequential.reserve(maxNbOfPairs);
  pairsOfNodesToProcess_all.reserve(maxNbOfPairs);

  // Search the pairs to merge and fill the vector
  for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
  {
    // Heuristic to determine with which adjacent node this current node has to merge.
    auto nodeIn = m_HeuristicFunc->GetBestAdjacentNode(&(*nodeIt));

    // The heuristic must return true if no adjacent node has been found.
    if(nodeIn != nullptr)
    {
      auto nodeOut = outputGraph->GetNodeAt(nodeIn->m_Edges.front().m_TargetId);

//      // Useless?
//      auto cost = nodeIn->m_Edges.front().m_Attributes.m_MergingCost;
//      (void) cost;

      // Both nodes must not have to be considered
      nodeIn->m_Valid = false;
      nodeOut->m_Valid = false;

      PairOfNodes newPair;
      newPair.first = nodeIn;
      newPair.second = nodeOut;

      // Keep trace of all pairs of nodes to merge
      pairsOfNodesToProcess_all.push_back(newPair);

      if (nodeIn->m_ThreadSafeForMerge && nodeOut->m_ThreadSafeForMerge)
        {
          // Add the pair of nodes to process in parallel
          pairsOfNodesToProcess_parallel.push_back(newPair);
        }
      else
        {
          // Add the pair of nodes to process in sequential
          pairsOfNodesToProcess_sequential.push_back(newPair);
        }

      merged = true;

    } // best adjacent node is not null

  } // next node

  // Store pre-merge processing time
  tmergePre.Stop(); timingsValues[3] += tmergePre.GetTotal();

  // Start parallel-merge processing time
  itk::TimeProbe tmergeMulti; tmergeMulti.Start();

  // Process the pairs of nodes in parallel
  int64_t nbOfPairs = pairsOfNodesToProcess_parallel.size();
  unsigned int nbOfThreads = this->GetNumberOfThreads();
  int64_t chunkSize = nbOfPairs/nbOfThreads;
  std::vector<std::thread> threadpoolMerge;
  threadpoolMerge.reserve(nbOfThreads);
  std::vector<std::pair<uint64_t,uint64_t>> ranges;
  for(unsigned int i=0; i < nbOfThreads; i++ )
    {
    // Compute range
    std::pair<uint64_t,uint64_t> range;
    range.first = i*chunkSize;
    if (i == nbOfThreads - 1 )
      {
      chunkSize += nbOfPairs % nbOfThreads;
      }
    range.second = range.first + chunkSize;
    ranges.push_back(range);

    // run thread
    threadpoolMerge.push_back(
        std::thread(FusionOfPairsInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
            std::ref( outputGraph ),
            std::ref( m_UpdateAttributeFunc ),
            std::ref( pairsOfNodesToProcess_parallel ),
            ranges[i].first,
            ranges[i].second
            )
    );

    } // next ProcessGraph() thread

  // Barrier for threadpoolMerge
  for (auto& t: threadpoolMerge) { t.join(); }

  // Stop parallel-merge processing time
  tmergeMulti.Stop(); timingsValues[4] += tmergeMulti.GetTotal();

  itk::TimeProbe tmergeSeq; tmergeSeq.Start();

  // merge non parallel
  FusionOfPairsInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>(
      outputGraph ,
      m_UpdateAttributeFunc,
      pairsOfNodesToProcess_sequential,
      0,
      pairsOfNodesToProcess_sequential.size()
  );

  tmergeSeq.Stop();

  // Stop merge processing time
  tmerge.Stop(); timingsValues[2] += tmerge.GetTotal();

  std::cout << tmerge.GetTotal() << "\t"
      << tmergePre.GetTotal() << "\t"
      << tmergeMulti.GetTotal() << "\t"
      << tmergeSeq.GetTotal() << "\t"
      << this->GetOutput()->GetNumberOfNodes() << std::endl;

  // Start remove() processing time
  itk::TimeProbe tremove; tremove.Start();

  outputGraph->RemoveNodes();

  // Stop remove() processing time
  tremove.Stop(); timingsValues[6] += tremove.GetTotal();

  // Store iteration processing time
  titeration.Stop();  timingsValues[0] += titeration.GetTotal();

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
	itk::TimeProbe recon_probe; recon_probe.Start();

  int64_t nbOfNodes = outputGraph->GetNumberOfNodes();
  unsigned int nbOfThreads = this->GetNumberOfThreads();
  int64_t chunkSize = nbOfNodes/nbOfThreads;

  std::vector<std::thread> threadpoolReconditioning;
  threadpoolReconditioning.reserve(nbOfThreads);
  std::vector<std::pair<uint64_t,uint64_t>> ranges;
  for(unsigned int i=0; i < nbOfThreads; i++ )
    {
    // Compute range
    std::pair<uint64_t,uint64_t> range;
    range.first = i*chunkSize;
    if (i == nbOfThreads - 1 )
      {
      chunkSize += nbOfNodes % nbOfThreads;
      }
    range.second = range.first + chunkSize;
    ranges.push_back(range);

    // run thread
    threadpoolReconditioning.push_back(
        std::thread(ReconditioningInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
            std::ref( outputGraph ),
            std::ref( m_MergingCostFunc ),
            ranges[i].first,
            ranges[i].second)
    );

    } // next ProcessGraph() thread

  // Barrier for threadpoolReconditioning
  for (auto& t: threadpoolReconditioning) { t.join(); }

	recon_probe.Stop();
	std::cout << "Reconditioning : " << recon_probe.GetTotal() << std::endl;

	// setup timings
	timingsValues.push_back(0.0f); timingsLabels.push_back("iterations");    // 0
  timingsValues.push_back(0.0f); timingsLabels.push_back("costs");         // 1
  timingsValues.push_back(0.0f); timingsLabels.push_back("merging");       // 2
  timingsValues.push_back(0.0f); timingsLabels.push_back("merging.pre");   // 3
  timingsValues.push_back(0.0f); timingsLabels.push_back("merging.multi"); // 4
  timingsValues.push_back(0.0f); timingsLabels.push_back("merging.post");  // 5
  timingsValues.push_back(0.0f); timingsLabels.push_back("remove");        // 6

	std::cout << "Iter.\t"
	    << "Costs\t"
      << "Merge\t"
      << "Merge(pre)\t"
      << "Merge(//)\t" // parallel
      << "Merge(--)\t" // sequential
	    << "nodes\t" << std::endl;

	uint64_t initialNbOfNodes = outputGraph->GetNumberOfNodes();
	this->UpdateProgress(.0f);
	for(uint32_t i = 0; i < m_MaxNumberOfIterations; i++)
	{
		std::cout << std::setprecision(5) << i+1 << "\t";
		if(!DoOneIteration())
		{
			m_MergingOver = true;
			break;
		}

		float progress = static_cast<float>(nbOfNodes);
		progress /= static_cast<float>(initialNbOfNodes);
		progress *= -1.0f;
		progress += 1.0f;
		this->UpdateProgress(progress);

	} // next iteration
	this->UpdateProgress(1.0f);
	std::cout << std::endl;

	// Display timings
	std::cout << " Overall processing time : " << std::endl;
	std::cout << "NbThreads";
  for (unsigned int i = 0 ; i < timingsLabels.size() ; i++)
    std::cout << "\t" <<  timingsLabels[i]  ;
  std::cout << std::endl;
  std::cout << this->GetNumberOfThreads();
  for (unsigned int i = 0 ; i < timingsLabels.size() ; i++)
    std::cout << "\t" << std::setprecision(5) << timingsValues[i];
  std::cout << std::endl;
  // Graph tming:

}



/*
 * Fully parallel
 */
template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc >
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::ComputeMergingCosts()
{

  // Retrieve the output graph.
  auto outputGraph = this->GetOutput();

  int64_t nbOfNodes = outputGraph->GetNumberOfNodes();
  unsigned int nbOfThreads = this->GetNumberOfThreads();
  int64_t chunkSize = nbOfNodes/nbOfThreads;

  std::vector<std::thread> threadpoolComputeMergingCosts, threadpoolReset;
  threadpoolComputeMergingCosts.reserve(nbOfThreads);
  threadpoolReset.reserve(nbOfThreads);
  std::vector<std::pair<uint64_t,uint64_t>> ranges;
  for(unsigned int i=0; i < nbOfThreads; i++ )
    {
    // Compute range
    std::pair<uint64_t,uint64_t> range;
    range.first = i*chunkSize;
    if (i == nbOfThreads - 1 )
      {
      chunkSize += nbOfNodes % nbOfThreads;
      }
    range.second = range.first + chunkSize;
    ranges.push_back(range);

    // run thread
    threadpoolComputeMergingCosts.push_back(
        std::thread(ComputeMergingCostsInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
            std::ref( outputGraph ),
            std::ref( m_MergingCostFunc ),
            ranges[i].first,
            ranges[i].second)
    );

    } // next ProcessGraph() thread

  // Barrier for threadpoolProcess
  for (auto& t: threadpoolComputeMergingCosts) { t.join(); }

  // Reset
  for(unsigned int i=0; i < nbOfThreads; i++ )
    {

    // run thread
    threadpoolReset.push_back(
        std::thread(ResetMergeFlagsInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
            std::ref( outputGraph ),
            std::ref( m_MergingCostFunc ),
            ranges[i].first,
            ranges[i].second)
    );

    } // next ResetGraph() thread

  // Barrier for threadpoolReset
  for (auto& t: threadpoolReset) { t.join(); }


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
