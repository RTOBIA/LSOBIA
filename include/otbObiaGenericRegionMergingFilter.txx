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
 *
 * TODO:
 * * something cleaner :]
 * * refactor the call to the worker, including the computation of chunkSize, ranges, ...
 * * maybe use itk thread classes (don't know if useful?)
 * * const stuff
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
  void SetPairsOfNodesToProcess(std::vector<std::pair<typename TOutputGraph::NodeType * ,
                                typename TOutputGraph::NodeType * >> &pairs)
  {
    m_PairsOfNodesToProcess = pairs;
  }
  void SetNumMergedNodes(std::vector<uint32_t> & vector)
  {
    m_NumMergedNodes = vector;
  }

  /*
   * Compute merging costs of the nodes
   * See "threadsafe.jpg" in the root directory for details.
   */
  void ComputeMergingCosts()
  {

    // Compute merging costs

    uint64_t minNodeId = m_OutputGraph->GetNumberOfNodes()+1, idx, minIdx;
    typename TMergingCostFunc::ValueType minCost;

    auto updateEdgesCosts = [this, &idx, &minIdx, &minCost, &minNodeId](typename TOutputGraph::EdgeType& edgeIt,
        typename TOutputGraph::NodeType& node){ // Lambda for the edges loop

      // Tell if this node has edges targeting nodes which are outside the thread range
      if ( edgeIt.m_TargetId < m_RangeStart || edgeIt.m_TargetId >= m_RangeEnd )
        {
        node.m_ThreadSafe = false;
        node.m_ThreadSafeForMerge = false;
        }

      // Retrieve the adjacent node.
      auto adjNode = m_OutputGraph->GetNodeAt(edgeIt.m_TargetId);
      if(m_MergingCostFunc->ComputeMergingCostsForThisAdjNode(adjNode))
        {

        // If one of the adjacent nodes
        // has merged at the previous iteration then we must compute the
        // merging cost.
        if(node.m_Attributes.m_HasPreviouslyMerged || adjNode->m_Attributes.m_HasPreviouslyMerged)
          {
          edgeIt.m_Attributes.m_MergingCost = m_MergingCostFunc->ComputeMergingCost(&node, adjNode);
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

    }; // end lambda for edges loop

    auto updateNodes = [this, &minNodeId, &minCost, &idx, &minIdx, &updateEdgesCosts](
        typename TOutputGraph::NodeType& node){ // Lambda for nodes loop
      node.m_ThreadSafe = true;
      node.m_ThreadSafeForMerge = true;

      if(m_MergingCostFunc->ComputeMergingCostsForThisNode(&(node))
          && node.m_Edges.size()>0) // Since the introducing of no-data, a node can be alone
        {

        // The merging cost function must give the maximum value of the merging cost.
        minCost = TMergingCostFunc::Max();
        idx = 0;
        minIdx = 0;
        node.m_HasToBeRemoved = false;
        node.m_Valid = true;

        // Loop over edges
        auto updateEdges = [&node, &updateEdgesCosts](typename TOutputGraph::EdgeType& edgeIt)  {
          updateEdgesCosts(edgeIt, node);
        };

        node.ApplyForEachEdge(updateEdges);

        // Finally we move the adjacent node with the lower merging cost
        // at the first position in the list of adjacent nodes.
        std::swap(node.m_Edges[0], node.m_Edges[minIdx]);
        }
    };

    m_OutputGraph->ApplyForEachNode(m_RangeStart, m_RangeEnd, updateNodes);

    // Grow the thread-unsafe-for-merge region (we know that the merge operation requires to R/W adjacent nodes of
    // the pair of nodes to merge)

    auto setThreadSafe = [this] (typename TOutputGraph::EdgeType& edge){
      if (m_RangeStart <= edge.m_TargetId && edge.m_TargetId < m_RangeEnd)
        {
        auto adjNode = m_OutputGraph->GetNodeAt(edge.m_TargetId);
        adjNode->m_ThreadSafeForMerge = false;
        }
    };

    auto growThreadRegions = [this, &setThreadSafe] (typename TOutputGraph::NodeType& node){
      if (node.m_ThreadSafe == false)
        {
        node.ApplyForEachEdge(setThreadSafe);
        }
    };

    m_OutputGraph->ApplyForEachNode(m_RangeStart, m_RangeEnd, growThreadRegions);

  } // ComputeMergingCosts()

  /*
   * Reset the merge flag(s) of nodes
   */
  void ResetMergeFlags()
  {
    auto setHasPrevMergedFalse = [] (typename TOutputGraph::NodeType& node){
      node.m_Attributes.m_HasPreviouslyMerged = false;
    };

    m_OutputGraph->ApplyForEachNode(m_RangeStart, m_RangeEnd, setHasPrevMergedFalse);
  } // ResetMergeFlags()

  /*
   * Performs the merge of pairs in the range [m_RangeStart, m_RangeEnd[
   */
  void FusionOfPairs()
  {
    for (auto pairIt =  m_PairsOfNodesToProcess.begin() + m_RangeStart ;
        pairIt != m_PairsOfNodesToProcess.begin() + m_RangeEnd; pairIt++)
      {
      typename TOutputGraph::NodeType * nodeIn = pairIt->first;
      typename TOutputGraph::NodeType * nodeOut = pairIt->second;

      // Update attributes of thread-safe nodes
      m_UpdateAttributesFunc->UpdateAttributes(nodeIn, nodeOut);

      // Merge the pair of nodes
      m_OutputGraph->MergePairOfNodes(nodeIn, nodeOut);
      }
  } // FusionOfPairs()

  /*
   * Reconditionning of the graph
   */
  void Reconditioning()
  {
    // Reset the edge attributes
    auto resetEdgAttr = [this](typename TOutputGraph::EdgeType& edge){
      edge.m_Attributes.m_MergingCost = m_MergingCostFunc->GetMax();
    };

    // Reset the node attributes
    auto resetNodeAttr = [&resetEdgAttr](typename TOutputGraph::NodeType& node){
      node.m_HasToBeRemoved = false;
      node.m_Valid = true;
      node.m_Attributes.m_HasPreviouslyMerged = true;
      node.ApplyForEachEdge(resetEdgAttr);
    };

    m_OutputGraph->ApplyForEachNode(m_RangeStart, m_RangeEnd, resetNodeAttr);
  }

  /*
   * Update nodes ids
   */
  void UpdateNodesIds()
  {
    // Update the edge target id
    auto lambdaDecrementIdEdge = [this](typename TOutputGraph::EdgeType& edge){
      edge.m_TargetId = edge.m_TargetId - m_NumMergedNodes[edge.m_TargetId];
    };

    // Update the node id
    auto lambdaDecrementIdNode = [this, &lambdaDecrementIdEdge](typename TOutputGraph::NodeType& node ){
      node.m_Id = node.m_Id - m_NumMergedNodes[node.m_Id];
      node.ApplyForEachEdge(lambdaDecrementIdEdge);
    };

    m_OutputGraph->ApplyForEachNode(m_RangeStart, m_RangeEnd, lambdaDecrementIdNode);
  }

private:
  TOutputGraph*           m_OutputGraph;
  TMergingCostFunc *      m_MergingCostFunc;
  TUpdateAttributeFunc *  m_UpdateAttributesFunc;
  int64_t                 m_RangeStart;
  int64_t                 m_RangeEnd;
  std::vector<std::pair<typename TOutputGraph::NodeType * ,typename TOutputGraph::NodeType * >>
  m_PairsOfNodesToProcess;
  std::vector<uint32_t>   m_NumMergedNodes;

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
void ResetMergeFlagsInRange (TOutputGraph * graph, uint64_t start, uint64_t end)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
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

template<typename TOutputGraph, typename TMergingCostFunc, typename TUpdateAttributeFunc>
void UpdateNodesInRange (TOutputGraph * graph, uint64_t start, uint64_t end,
                         std::vector<uint32_t> & numMergedNodes)
{

  ThreadWorker<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc> worker;
  worker.SetGraph( graph );
  worker.SetRange(start, end);
  worker.SetNumMergedNodes(numMergedNodes);
  worker.UpdateNodesIds();

}

template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc>
bool
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::DoOneIteration()
{
  // Output graph
  auto outputGraph = this->GetOutput();

  // Measure iteration total time
  itk::TimeProbe titeration;  titeration.Start();

  // Compute the merging costs for all the pairs of adjacent nodes.
  itk::TimeProbe tcosts;  tcosts.Start();
  ComputeMergingCosts();

  tcosts.Stop(); timingsValues[1] += tcosts.GetTotal();
  std::cout << std::setprecision(3) << tcosts.GetTotal() << "\t";

  // Flag indicating if at least one merge has been done during the iteration.
  bool merged = false;

  itk::TimeProbe tmerge;  tmerge.Start();
  itk::TimeProbe tmergePre;  tmergePre.Start();

  // First part is not multithreaded, we just identify pairs of nodes to merge
  // TODO: this 3 std::vector of pairs of nodes pointers seems not the best option !
  //       pairsOfNodesToProcess_all has to be kept, but maybe the following vectors
  //       might just be vectors or pointers:
  //        * pairsOfNodesToProcess_parallel
  //        * pairsOfNodesToProcess_sequential
  typedef std::pair<NodeType * ,NodeType * > PairOfNodes;
  std::vector<PairOfNodes> pairsOfNodesToProcess_parallel;
  std::vector<PairOfNodes> pairsOfNodesToProcess_sequential;
  std::vector<PairOfNodes> pairsOfNodesToProcess_all;

  // Let's allocate the maximum number of pairs possible (ie 1+n/2)
  const int64_t maxNbOfPairs = outputGraph->GetNumberOfNodes() / 2 + 1;
  pairsOfNodesToProcess_parallel.reserve(maxNbOfPairs);
  pairsOfNodesToProcess_sequential.reserve(maxNbOfPairs);
  pairsOfNodesToProcess_all.reserve(maxNbOfPairs);

  // Lambda to search the pairs to merge and fill the vector
  auto searchPairs = [this, &outputGraph, &pairsOfNodesToProcess_parallel, &pairsOfNodesToProcess_sequential,
                      &pairsOfNodesToProcess_all, &merged] (typename TOutputGraph::NodeType& node ){

    // Heuristic to determine with which adjacent node this current node has to merge.
    auto nodeIn = m_HeuristicFunc->GetBestAdjacentNode(&(node));

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

      // TODO: push_back() seems to copy the input.
      //       Maybe something smarted is required here (vector of pointers?).
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

  };

  outputGraph->ApplyForEachNode(searchPairs);

  // Store pre-merge processing time
  tmergePre.Stop(); timingsValues[3] += tmergePre.GetTotal();

  // Start parallel-merge processing time
  itk::TimeProbe tmergeMulti; tmergeMulti.Start();

  // Process the pairs of nodes in parallel
  const int64_t nbOfPairs = pairsOfNodesToProcess_parallel.size();
  const unsigned int nbOfThreads = this->GetNumberOfThreads();
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

  //  // Remove nodes
  //  outputGraph->RemoveNodes();

  //////////////////////////////////////  
  // Remove nodes (in parallel)
  //////////////////////////////////////
  std::vector<uint32_t> numMergedNodes = outputGraph->RemoveNodes(false);

  const int64_t nbOfNodes = outputGraph->GetNumberOfNodes();
  chunkSize = nbOfNodes/nbOfThreads;
  std::vector<std::thread> threadpoolUpdateNodes;
  threadpoolUpdateNodes.reserve(nbOfThreads);
  ranges.clear();
  for(unsigned int i = 0; i < nbOfThreads; i++ )
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

    // run threads to count nodes to remove
    threadpoolUpdateNodes.push_back(
        std::thread(UpdateNodesInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
                    std::ref( outputGraph ),
                    ranges[i].first,
                    ranges[i].second,
                    std::ref( numMergedNodes))
    );

    } // next ProcessGraph() thread

  // Barrier for threadpoolProcess
  for (auto& t: threadpoolUpdateNodes) { t.join(); }

  //////////////////////////////////////

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

  const int64_t nbOfNodes = outputGraph->GetNumberOfNodes();
  const unsigned int nbOfThreads = this->GetNumberOfThreads();
  int64_t chunkSize = nbOfNodes/nbOfThreads;

  std::vector<std::thread> threadpoolReconditioning;
  threadpoolReconditioning.reserve(nbOfThreads);
  std::vector<std::pair<uint64_t,uint64_t>> ranges;
  for(unsigned int i = 0; i < nbOfThreads; i++ )
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

  this->UpdateProgress(.0f);
  for(uint32_t i = 0; i < m_MaxNumberOfIterations; i++)
    {
    std::cout << std::setprecision(5) << i+1 << "\t";
    if(!DoOneIteration())
      {
      m_MergingOver = true;
      break;
      }

    // TODO: it should be possible to estimate the progress
    // with a model (a,b) like n(k)=exp(-k.a)+b where n(k) is the number
    // of nodes at each iteration k
    //		float progress = ...
    //		this->UpdateProgress(progress);

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

}

/*
 * Compute merging costs.
 * This function is fully parallelized.
 */
template< typename TInputGraph, typename TOutputGraph, typename TMergingCostFunc, typename THeuristic, typename TUpdateAttributeFunc >
void
GenericRegionMergingFilter<TInputGraph, TOutputGraph, TMergingCostFunc, THeuristic, TUpdateAttributeFunc>::ComputeMergingCosts()
{

  // Retrieve the output graph.
  auto outputGraph = this->GetOutput();

  const int64_t nbOfNodes = outputGraph->GetNumberOfNodes();
  const unsigned int nbOfThreads = this->GetNumberOfThreads();
  int64_t chunkSize = nbOfNodes/nbOfThreads;

  std::vector<std::thread> threadpoolComputeMergingCosts, threadpoolReset;
  threadpoolComputeMergingCosts.reserve(nbOfThreads);
  threadpoolReset.reserve(nbOfThreads);
  std::vector<std::pair<uint64_t,uint64_t>> ranges;
  for(unsigned int i = 0; i < nbOfThreads; i++ )
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
  for(unsigned int i = 0; i < nbOfThreads; i++ )
    {
    // run thread
    threadpoolReset.push_back(
        std::thread(ResetMergeFlagsInRange<TOutputGraph, TMergingCostFunc, TUpdateAttributeFunc>,
                    std::ref( outputGraph ),
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
