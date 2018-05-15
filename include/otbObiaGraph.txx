#ifndef otbObiaGraph_txx
#define otbObiaGraph_txx
#include "otbObiaGraph.h"
namespace otb
{

namespace obia
{

template< typename TNode >
uint64_t
Edge<TNode>::GetNumberOfBytesToSerialize() const
{
    return CoordValueSize + // The starting coordinates of the target node
           UInt32Size + // the boundary length
           m_Attributes.GetNumberOfBytesToSerialize();
}

template< typename TNodeAttribute, typename TEdgeAttribute >
typename Node<TNodeAttribute, TEdgeAttribute>::EdgeType* 
Node<TNodeAttribute, TEdgeAttribute>::AddEdge()
{
    m_Edges.emplace_back(); 
    return &(m_Edges.back());
}

template< typename TNodeAttribute, typename TEdgeAttribute >
typename Node<TNodeAttribute, TEdgeAttribute>::EdgeIteratorType
Node<TNodeAttribute, TEdgeAttribute>::FindEdge(const IdType targetId)
{
    return std::find_if(m_Edges.begin(), m_Edges.end(), [&](const EdgeType& edge)->bool{
            return edge.m_TargetId == targetId;
    });
}

template< typename TNodeAttribute, typename TEdgeAttribute >
uint64_t 
Node<TNodeAttribute, TEdgeAttribute>::GetMemorySize() const
{
  uint64_t nodeMemory = sizeof(Self);

  nodeMemory+=m_Contour.GetMemorySize()+m_Attributes.GetMemorySize();

  for(auto edgeIt = m_Edges.begin(); edgeIt != m_Edges.end(); edgeIt++)
    {
        nodeMemory += edgeIt->GetMemorySize();
    }


    return nodeMemory;
}

template< typename TNodeAttribute, typename TEdgeAttribute >
uint64_t
Node<TNodeAttribute, TEdgeAttribute>::GetNumberOfBytesToSerialize() const
{
    uint64_t numberOfBytes = 0;

    // Bounding box number of bytes
    numberOfBytes += m_BoundingBox.size() * UInt32Size;

    // Contour number of bytes
    numberOfBytes += m_Contour.GetNumberOfBytesToSerialize();

    // Specific attributes number of bytes
    numberOfBytes += m_Attributes.GetNumberOfBytesToSerialize();

    // Number of edges
    numberOfBytes += UInt32Size;

    for(auto& edg : m_Edges)
    {
        numberOfBytes += edg.GetNumberOfBytesToSerialize();
    }

    return numberOfBytes;
}

template< typename TNode >
void 
Graph<TNode>::RemoveEdgesToThisNode(NodeType& node)
{
    // Goal: remove the edges targeting to this node
    auto lambdaOp = [&](EdgeType& edg){

        // Retrieve the adjacent node of this node
        auto adjNode = GetNodeAt(edg.m_TargetId);

        // Find the outgoing edge from adjNode to node
        // This node has to exist !
        auto adjToNode = adjNode->FindEdge(node.m_Id);

        // We can safely remove this edge
        adjNode->m_Edges.erase(adjToNode);
    };

    node.ApplyForEachEdge(lambdaOp);    

    // The node can be marked as to be removed.
    node.m_HasToBeRemoved = true;
}

template< typename TNode >
void 
Graph<TNode>::SetNumberOfNodes(const uint64_t numNodes)
{
    m_Nodes.reserve(numNodes);
}

template< typename TNode >
typename Graph<TNode>::NodeType* 
Graph<TNode>::AddNode()
{
    m_Nodes.emplace_back();
    return &(m_Nodes.back());
}

template< typename TNode >
void
Graph<TNode>::InitStartingNode(NodeType* node, const IdType id)
{
    // Initialisation of its position
    node->m_Id = id;
        
    // Initialisation of the bounding box
    node->m_BoundingBox[0] = id % m_ImageWidth;
    node->m_BoundingBox[1] = id / m_ImageWidth;
    node->m_BoundingBox[2] = 1;
    node->m_BoundingBox[3] = 1;

    // Initialisation of the contour
    node->m_Contour.SetStartingCoords(id);
    node->m_Contour.FirstInit();

    // We don't fo this here anymore.
    // Now we use nodata to init properly the edges in the imageToGraph
    
//    // Initialisation of the edges
//    auto neighbors = otb::obia::SpatialTools::FourConnectivity(id, m_ImageWidth, m_ImageHeight);
//
//    uint32_t numEdges = 0;
//    for(unsigned short n = 0; n < 4; n++){ if(neighbors[n] > -1){ numEdges++; } }
//    node->m_Edges.reserve(numEdges);
//
//    for(unsigned short n = 0; n < 4; n++)
//    {
//        if(neighbors[n] > -1)
//        {
//            // Add an edge to the current node targeting the adjacent node
//            auto newEdge = node->AddEdge();
//
//            // Add the target
//            newEdge->m_TargetId = neighbors[n];
//
//            // Initialisation of the boundary
//            newEdge->m_Boundary = 1;
//        }
//    }
}

template< typename TNode >
uint64_t 
Graph<TNode>::GetNumberOfNodes() const
{
    return m_Nodes.size();
}

template< typename TNode >
typename Graph<TNode>::NodeType* 
Graph<TNode>::GetNodeAt(const IdType id)
{
	if(id > (m_Nodes.size() - 1) || id < 0)
	{
		std::cout <<"NULL NODE FOR " << id << std::endl;
		exit(EXIT_FAILURE);
		return nullptr;
	}
    return &(m_Nodes[id]);
}

template< typename TNode >
void 
Graph<TNode>::Print()
{
    for(auto i = 0; i < GetNumberOfNodes(); i++)
    {
        auto node = GetNodeAt(i);
        std::cout << node->GetFirstPixelCoords() << "(" << node->m_Id << "): ";
        for(auto edgeIt = node->m_Edges.begin(); edgeIt != node->m_Edges.end(); edgeIt++)
        {
            std::cout << (GetNodeAt(edgeIt->m_TargetId))->GetFirstPixelCoords() << "("<< edgeIt->m_TargetId << ") \t";
        }
        std::cout << std::endl;
    }
}

template< typename TNode >
void 
Graph<TNode>::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

template< typename TNode >
template< typename LambdaFunctionType >
void Graph<TNode>::ApplyForEachNode(LambdaFunctionType f)
{
    std::for_each(m_Nodes.begin(), m_Nodes.end(), f);
}

template< typename TNode >
template< typename LambdaFunctionType >
void Graph<TNode>::ApplyForEachNode(const uint64_t rangeStart, const uint64_t rangeEnd, LambdaFunctionType f)
{
    std::for_each(m_Nodes.begin()+rangeStart, m_Nodes.begin()+rangeEnd, f);
}

template< typename TNode >
void
Graph<TNode>::MergeEdge(NodeType* nodeIn, NodeType* nodeOut)
{

  // Explore the edges of nodeOut
  for(auto& edgeIt : nodeOut->m_Edges)
    {
    // Local variable to record the boundary length of nodeOut with its
    // adjacent node.
    uint32_t boundary;

    // Retrieve the adjacent node of nodeOut -> adjNodeOut
    auto adjNodeOut = GetNodeAt(edgeIt.m_TargetId);

    // Find the edge from adjNodeOut to nodeOut
    auto edgAdjToOut = adjNodeOut->FindEdge(nodeOut->m_Id);

    // If this edge is the first edge of adjNodeOut then adjNodeOut is not valid
    // anymore for merging with one of its adjacent node during a segmentation process.
    // This is not so generic since some obia processes won't need that but it is a price
    // to pay to factorize and shorten a bit the code.
    if(edgAdjToOut == adjNodeOut->m_Edges.begin())
      {
      adjNodeOut->m_Valid = false;
      }

    // Record the boundary length between adjNodeOut and nodeOut
    boundary = edgAdjToOut->m_Boundary;

    // The edge adjNodeOut -> nodeOut can be safely removed
    adjNodeOut->m_Edges.erase(edgAdjToOut);

    // Edges have to be added or updated if adjNodeOut is not nodeIn
    if(adjNodeOut != nodeIn)
      {
      // Try to find if there is an edge between adjNodeOut and nodeIn
      auto edgAdjToIn = adjNodeOut->FindEdge(nodeIn->m_Id);

      if(edgAdjToIn == adjNodeOut->m_Edges.end())
        {
        // If the edge does not exist, it has to be created

        // AdjNodeOut -> nodeIn
        auto adjOutToIn = adjNodeOut->AddEdge();
        adjOutToIn->m_TargetId = nodeIn->m_Id;
        adjOutToIn->m_Boundary = boundary;

        // nodeIn -> AdjNodeOut
        auto inToAdjOut = nodeIn->AddEdge();
        inToAdjOut->m_TargetId = adjNodeOut->m_Id;
        inToAdjOut->m_Boundary = boundary;

        } // end if(edgAdjToIn == adjNodeOut->m_Edges.end())
      else
        {
        // the edges exist, the boundary attribute has to be updated

        // Increment AdjNodeOut -> nodeIn
        edgAdjToIn->m_Boundary += boundary;

        // Increment nodeIn -> AdjNodeOut
        auto edgInToAdj = nodeIn->FindEdge(adjNodeOut->m_Id);
        edgInToAdj->m_Boundary += boundary;
        }

      } // end if(adjNodeOut != nodeIn)

    } // end for(auto edgeIt = nodeOut->m_Edges.begin(); edgeIt != nodeOut->m_Edges.end(); edgeIt++)

  // All the edges of nodeOut can be safely removes
  typename NodeType::EdgeListType().swap(nodeOut->m_Edges);
}

template< typename TNode >
void
Graph<TNode>::MergePairOfNodes(NodeType* nodeIn, NodeType* nodeOut)
{
  // Fusion of the bounding box
  SpatialTools::MergeBoundingBox(nodeIn->m_BoundingBox, nodeOut->m_BoundingBox);

  // Fusion of the contour
  nodeIn->m_Contour.MergeWith(nodeOut->m_Contour, m_ImageWidth, m_ImageHeight);

  // Fusion of the edges
  MergeEdge(nodeIn, nodeOut);

  // nodeOut has to be removed, so we mark it as it
  nodeOut->m_HasToBeRemoved = true;
}

template< typename TNode >
void 
Graph<TNode>::Merge(NodeType* nodeIn, NodeType* nodeOut)
{
  // Merge the pair of nodes
  MergePairOfNodes(nodeIn, nodeOut);

  // Both nodes must not have to be considered
  nodeIn->m_Valid = false;
  nodeOut->m_Valid = false;
}

template< typename TNode >
std::vector<uint32_t> 
Graph<TNode>::RemoveNodes(bool update)
{
  // To keep nodes aligned in memory, we avoid using vector of pointers
  // but directly vector of nodes. However this complicates the removal
  // operation since the ids of the nodes has to be updated.
  // Space complexity: O(Nt + Et) where Nt is the total number of nodes and Et the total number of edges
  // Time complexity: O(Nt * E) where Nt is the total number of nodes and E the average number of edges per node.

  // Count the number of merged nodes at each position
  // Linear in time wrt to the number of nodes
  std::vector<uint32_t> numMergedNodes(m_Nodes.size(), 0);
  uint64_t idx = 0;
  uint32_t numMerged = 0;

  auto lambdaNumMerges = [&numMergedNodes, &idx, &numMerged](const NodeType& node){

    if(node.m_HasToBeRemoved)
      {
      numMerged++;
      }

    numMergedNodes[idx] = numMerged;
    idx++;
  };

  ApplyForEachNode(lambdaNumMerges);

  // Remove the nodes: Linear in time wrt to the number of nodes
  auto eraseIt = std::remove_if(m_Nodes.begin(), m_Nodes.end(), [](NodeType& node){
    return node.m_HasToBeRemoved == true;
  });

  m_Nodes.erase(eraseIt, m_Nodes.end());

  // We might want do this in parallel instead
  if (update)
    {
    auto lambdaDecrementIdEdge = [&numMergedNodes](EdgeType& edge){
      edge.m_TargetId = edge.m_TargetId - numMergedNodes[edge.m_TargetId];
    };

    auto lambdaDecrementIdNode = [&numMergedNodes, &lambdaDecrementIdEdge]( NodeType& node ){
      // Update the node's id
      node.m_Id = node.m_Id - numMergedNodes[node.m_Id];
      node.ApplyForEachEdge(lambdaDecrementIdEdge);
    };

    ApplyForEachNode(lambdaDecrementIdNode);
    }   
  return numMergedNodes;

}

template< typename TNode >
uint64_t
Graph<TNode>::GetMemorySize() const
{
	unsigned int nbImageInfos = 5;
    uint64_t graphMemory = nbImageInfos * UInt32Size + // image infos
                           sizeof(std::string) + // container for projection
                           m_ProjectionRef.size() * CharSize +
                           sizeof(NodeListType);

    for(auto& node : m_Nodes)
    {
        graphMemory += node.GetMemorySize();
    }

    return graphMemory;
}


} // end of namespace obia

} // end of namespace otb

#endif
