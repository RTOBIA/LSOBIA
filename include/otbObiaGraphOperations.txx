#ifndef otbObiaGraphOperations_txx
#define otbObiaGraphOperations_txx
#include "otbObiaGraphOperations.h"
#include "otbObiaStreamUtils.h"

//TEMP
#include "otbMPIConfig.h"

namespace otb
{
namespace obia
{

template< typename TGraph >
void 
GraphOperations<TGraph>::RemoveUnstableNodes(GraphPointerType graph,
                                             const ProcessingTile& tile,
                                             const uint32_t imageWidth)
{
    // Lower and upper bounds of the tile without the margin
    const uint32_t lowerCol = tile.m_Frame.GetIndex(0);
    const uint32_t lowerRow = tile.m_Frame.GetIndex(1);
    const uint32_t upperCol = lowerCol + tile.m_Frame.GetSize(0) - 1;
    const uint32_t upperRow = lowerRow + tile.m_Frame.GetSize(1) - 1;

    // Flag indicating if the node has to be kept in the graph.
    bool stableNode;

    // Coordinates of a pixel
    uint32_t x, y;

    auto lambdaOp = [&](NodeType& node){

        // First assumption: the node is stable.
        stableNode = true;

        if(SpatialTools::IsBboxOutsideBoundaries(node.m_BoundingBox,
                                                   lowerCol,
                                                   lowerRow,
                                                   upperCol,
                                                   upperRow))
        {
            // The node is outside, it is obviously unstable.
            stableNode = false;
        }
        else if(!SpatialTools::IsBboxInsideBoundaries(node.m_BoundingBox,
                                                       lowerCol,
                                                       lowerRow,
                                                       upperCol,
                                                       upperRow))
        {
            // Second assumption: the node is not stable
            stableNode = false;

            // We generate the border pixels located at the exterior ring of the node
            std::unordered_set<CoordValueType> borderPixels;
            node.m_Contour.GenerateBorderPixels(borderPixels, imageWidth);

            // Loop over the pixels
            for(auto& pix : borderPixels)
            {
                x = pix % imageWidth;
                y = pix / imageWidth;

                if( x >= lowerCol && x <= upperCol && y >= lowerRow && y <= upperRow)
                {
                    // If the node has at least one pixel in the stable area
                    // then the node is stable.
                    stableNode = true;

                    // The loop can be stopped.
                    break;
                } // end if( x >= lowerCol && x <= upperCol && y >= lowerRow && y <= upperRow)

            } // end for(auto& pix : borderPixels)

        } // end else if

        if(stableNode == false)
        {
            graph->RemoveEdgesToThisNode(node);
        }

    }; // end lambdaOp

    // Apply the lambda function to each node.
    graph->ApplyForEachNode(lambdaOp);

    // Remove all unstable nodes from the graph
    graph->RemoveNodes();

}

template< typename TGraph >
typename GraphOperations<TGraph>::GraphPointerType
GraphOperations<TGraph>::DeSerializeGraph(const std::vector<char>& serializedGraph)
{
    // Local variable: offset in the bit stream.
    uint64_t position = 0;

    // Initialization of the resulting graph
    auto graph = GraphType::New();

    // Read the number of nodes
    uint64_t numNodes;
    from_stream(serializedGraph,numNodes,position);
    

    // Reserve memory space for the graph
    graph->SetNumberOfNodes(numNodes);

    // Map: starting starting coords of each node -> its id
    std::unordered_map< CoordValueType, IdType > startPixIdMap;

    // Loop over the nodes
    for(IdType id = 0; id < numNodes; id++)
    {

        // Add a new node at the end of the adjacency list
        // of the graph.
        auto newNode = graph->AddNode();

        // The id of the node
        newNode->m_Id = id;

        // Set the flags appropriately
        newNode->m_HasToBeRemoved = false;
        newNode->m_Valid = true;

        // The bounding box
	from_stream(serializedGraph,newNode->m_BoundingBox,position);

        // The contour
        newNode->m_Contour.DeSerialize(serializedGraph, position);

        // Create a new entry in the map startPixIdMap
        assert(startPixIdMap.find(newNode->GetFirstPixelCoords()) == startPixIdMap.end());

	startPixIdMap[newNode->GetFirstPixelCoords()] = id;
	
        // The specific attributes of the node
        newNode->m_Attributes.DeSerialize(serializedGraph, position);

        // The number of edges
        uint32_t numEdges;
	from_stream(serializedGraph,numEdges,position);

        newNode->m_Edges.reserve(numEdges);

        for(uint32_t e = 0; e < numEdges; e++)
        {
            // Add a new edge
            auto newEdge = newNode->AddEdge();

            // The starting coords
	    from_stream(serializedGraph,newEdge->m_TargetId,position);
	    
	    // The boundary
	    from_stream(serializedGraph,newEdge->m_Boundary,position);

            // Serialize the specific attributes of the edges
            newEdge->m_Attributes.DeSerialize(serializedGraph, position);

        } // end for(uint32_t e = 0; e < numEdges; e++) 

    } // end for(IdType id = 0; id < numNodes; id++)

    // Update the target ids of the edges with the real ids of the nodes
    for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    {
        for(auto& edg : nodeIt->m_Edges)
        {
            auto isInStartPixIdMap = startPixIdMap.find(edg.m_TargetId);

	    assert(isInStartPixIdMap != startPixIdMap.end());
	    edg.m_TargetId = isInStartPixIdMap->second;

        } // end for(auto& edg : nodeIt->m_Edges)

    } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    
    return graph;
}

template< typename TGraph >
std::vector< char >
GraphOperations<TGraph>::SerializeStabilityMargin(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                                                    const GraphPointerType graph)
{
    // The resulting serialized stability margin of the graph: it is jus a sequence of bytes.
    std::vector< char > serializedGraph;

    // The first pass consists of determining the number of bytes to write.
    std::vector< uint32_t > numEdgesPerNode(stabilityMargin.size(), 0);

    // Number of nodes
    uint64_t numberOfBytes = UInt64Size;

    // Running id.
    IdType i = 0;

    // Loop over the map: node pointer -> exploration depth (useless here)
    for(const auto& kv : stabilityMargin)
    {
        // Retrieve the node.
        auto node = kv.first;

        // Serialization of the bounding box
        numberOfBytes += stream_offset(node->m_BoundingBox);

        // Serialization of the contour
        numberOfBytes += node->m_Contour.GetNumberOfBytesToSerialize();

        // Serialization of the specific attributes
        numberOfBytes += node->m_Attributes.GetNumberOfBytesToSerialize();

        // Serialization of the number of edges
        numberOfBytes += stream_offset(node->m_Edges.size());

        // Loop over the edges
        for(const auto& edg : node->m_Edges)
        {

            // Is it a dangling edge ? i.e an edge that targets a node which
            // is not in the stability margin.
            if(stabilityMargin.find(graph->GetNodeAt(edg.m_TargetId)) != stabilityMargin.end())
            {
                  numEdgesPerNode[i]++;
                  numberOfBytes += edg.GetNumberOfBytesToSerialize();
            }

        } // end for(const auto& edg : node->m_Edges)

        // increment the id
        i++;

    } // end for(const auto& kv : stabilityMargin)

    // Second pass: build the serialized graph
    serializedGraph.assign(numberOfBytes, char());

    // Reset the id value
    i = 0;

    // Displacement value in the sequence of bytes.
    uint64_t position = 0;

    // Serialize the number of nodes.
    const uint64_t numNodes = stabilityMargin.size();
    to_stream(serializedGraph,numNodes,position);

    for(const auto& kv : stabilityMargin)
    {
        auto node = kv.first;

        // Serialization of the bounding box
	to_stream(serializedGraph,node->m_BoundingBox,position);

        // Serialization of the contour
        node->m_Contour.Serialize(serializedGraph, position);

        // Serialization of the specific attributes
        node->m_Attributes.Serialize(serializedGraph, position);

        // Serialization of the number of edges
	to_stream(serializedGraph,numEdgesPerNode[i],position);

        for(const auto& edg : node->m_Edges)
        {
              if(stabilityMargin.find(graph->GetNodeAt(edg.m_TargetId)) != stabilityMargin.end())
            {
              
              // Serialize the first starting coords
              const CoordValueType startPix = (graph->GetNodeAt(edg.m_TargetId))->GetFirstPixelCoords();
	      to_stream(serializedGraph,startPix,position);

              // Serialize the boundary
	      to_stream(serializedGraph,edg.m_Boundary,position);

              // Serialize the specific attributes of the edges
              edg.m_Attributes.Serialize(serializedGraph, position);

            } // end if(stabilityMargin.find(graph->GetNodeAt(edg.m_TargetId)) != stabilityMargin.end())

        } // end for(const auto& edg : node->m_Edges)

        i++;

    } // end for(const auto& kv : stabilityMargin)


  return serializedGraph;
}

template< typename TGraph >
std::vector<char>
GraphOperations<TGraph>::SerializeGraph(const GraphPointerType graph)
{
    std::vector< char > serializedGraph;

    // The first pass consists of determining the number of bytes to write.

    // Number of nodes
    uint64_t numberOfBytes = stream_offset(graph->GetNumberOfNodes());

    for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    {
        // Serialization of the bounding box
        numberOfBytes += stream_offset(nodeIt->m_BoundingBox);

        // Serialization of the contour
        numberOfBytes += nodeIt->m_Contour.GetNumberOfBytesToSerialize();

        // Serialization of the specific attributes
        numberOfBytes += nodeIt->m_Attributes.GetNumberOfBytesToSerialize();

        // Serialization of the number of edges
	numberOfBytes+= stream_offset(nodeIt->m_Edges.size());

        for(const auto& edg : nodeIt->m_Edges)
        {

          numberOfBytes += edg.GetNumberOfBytesToSerialize();

        } // end for(const auto& edg : nodeIt->m_Edges)

    } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)

    serializedGraph.assign(numberOfBytes, char());

    // Second pass: build the serialized graph
    uint64_t position = 0;
    const uint64_t numNodes = graph->GetNumberOfNodes();

    to_stream(serializedGraph,numNodes,position);

    for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    {
        // Serialization of the bounding box
        to_stream(serializedGraph,nodeIt->m_BoundingBox,position);
        
        // Serialization of the contour
        nodeIt->m_Contour.Serialize(serializedGraph, position);

        // Serialization of the specific attributes
        nodeIt->m_Attributes.Serialize(serializedGraph, position);

        // Serialization of the number of edges
        const uint32_t numEdges = nodeIt->m_Edges.size();
	to_stream(serializedGraph,numEdges,position);

        for(const auto& edg : nodeIt->m_Edges)
        {
            // Serialize the first starting coords
             const CoordValueType startPix = (graph->GetNodeAt(edg.m_TargetId))->GetFirstPixelCoords();
	     to_stream(serializedGraph,startPix,position);
	    //to_stream(serializedGraph,edg.m_TargetId,position);

            // Serialize the boundary
	    to_stream(serializedGraph,edg.m_Boundary,position);

            // Serialize the specific attributes of the edges
            edg.m_Attributes.Serialize(serializedGraph, position);

        } // end for(const auto& edg : nodeIt->m_Edges)

    } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)

    return serializedGraph;
}

template< typename TGraph >
void
GraphOperations<TGraph>::WriteMarginGraphToDisk(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                        const GraphPointerType graph,
                        const std::string& outputPath)
{
    // Serialize the stability margin
    const auto serializedStabilityMargin = SerializeStabilityMargin(stabilityMargin, graph);

    // Open the file stream
    std::ofstream outFile(outputPath, std::ios::out | std::ios::binary);
    
    if(outFile.good())
    {
        // Retrieve the number of bytes to write.
        uint64_t numberOfBytes = serializedStabilityMargin.size();
        outFile.write(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);
        outFile.write(&serializedStabilityMargin[0], serializedStabilityMargin.size() * CharSize);
        outFile.close();
    }
    else
    {
        std::cerr << "Cannot open ouput file to write the graph: " << outputPath << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename TGraph >
void
GraphOperations<TGraph>::WriteSerializedMarginToDisk(const std::vector<char>& serializedStabilityMargin,
                         const std::string& outputPath)
{
    // Open the file stream
    std::ofstream outFile(outputPath, std::ios::out | std::ios::binary);

    if(outFile.good())
    {
        // Retrieve the number of bytes to write.
        uint64_t numberOfBytes = serializedStabilityMargin.size();
        outFile.write(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);
        outFile.write(&serializedStabilityMargin[0], serializedStabilityMargin.size() * CharSize);
        outFile.close();
    }
    else
    {
        std::cerr << "Cannot open ouput file to write the graph: " << outputPath << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename TGraph >
std::vector< char >
GraphOperations<TGraph>::ReadSerializedMarginFromDisk(const std::string& inputPath)
{
    // Open the file stream
    std::ifstream inFile(inputPath, std::ios::in | std::ios::binary);

    if(inFile.good())
    {
        // Retrieve the number of bytes to read.
        uint64_t numberOfBytes;
        inFile.read(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);
        std::vector<char> serializedGraph(numberOfBytes);
        inFile.read(&serializedGraph[0], numberOfBytes * CharSize);
        inFile.close();
        return serializedGraph;
    }
    else
    {
        std::cerr << "Cannot open input file to load the graph: " << inputPath << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename TGraph >
void
GraphOperations<TGraph>::WriteGraphToDisk(const GraphPointerType graph,
                  const std::string& outputPath)
{
    // Serialize the graph
    const auto serializedGraph = SerializeGraph(graph);

    // Open the file stream
    std::ofstream outFile(outputPath, std::ios::out | std::ios::binary  | std::ios::trunc);

    if(outFile.good())
    {
        // Retrieve the number of bytes to write.
        uint64_t numberOfBytes = serializedGraph.size();
        outFile.write(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);
        outFile.write(&serializedGraph[0], serializedGraph.size() * CharSize);
        outFile.close();
    }
    else
    {
        std::cerr << "Cannot open ouput file to write the graph: " << outputPath << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Ecriture du header */
    GraphOperations<TGraph>::WriteGraphHeader(outputPath, graph);
}

template< typename TGraph >
typename GraphOperations<TGraph>::GraphPointerType
GraphOperations<TGraph>::ReadGraphFromDisk(const std::string & inputPath)
{
    // Open the file stream
    std::ifstream inFile(inputPath, std::ios::in | std::ios::binary);

    if(inFile.good())
    {
        // Retrieve the number of bytes to read.
        uint64_t numberOfBytes;
        inFile.read(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);
        std::vector<char> serializedGraph(numberOfBytes);
        inFile.read(&serializedGraph[0], numberOfBytes * CharSize);
        inFile.close();
        
        auto graph = DeSerializeGraph(serializedGraph);

        /* Lecture du header et update du graphe avec les informations */
        GraphImageInfo info = GraphOperations<TGraph>::ReadGraphHeader(inputPath);
        graph->SetImageWidth(info.m_width);
        graph->SetImageHeight(info.m_height);
        graph->SetNumberOfSpectralBands(info.m_nbBands);
        graph->SetOriginX(info.m_originX);
        graph->SetOriginY(info.m_originY);
        graph->SetProjectionRef(info.m_projectionRef);
        return graph;
    }
    else
    {
        std::cerr << "Cannot open input file to load the graph: " << inputPath << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename TGraph >
std::vector< typename GraphOperations<TGraph>::NodeType* >
GraphOperations<TGraph>::GetListOfBorderNodes(const GraphPointerType graph,
                                              const ProcessingTile& tile,
                                              const uint32_t nbTilesX,
                                              const uint32_t nbTilesY,
                                              const uint32_t inputLSImageWidth)
{
    // Lower and upper bounds of the tile without the margin
    // Lower and upper bounds of the tile without the margin
    const uint32_t lowerCol = tile.m_Frame.GetIndex(0);
    const uint32_t lowerRow = tile.m_Frame.GetIndex(1);
    const uint32_t upperCol = lowerCol + tile.m_Frame.GetSize(0) - 1;
    const uint32_t upperRow = lowerRow + tile.m_Frame.GetSize(1) - 1;

    // This list contains the border nodes.
    std::vector<NodeType*> listOfBorderNodes;

    auto lambdaOp = [&](NodeType& node){

	if(!SpatialTools::IsBboxStrictlyInsideBoundaries(node.m_BoundingBox,
													 lowerCol,
													 lowerRow,
													 upperCol,
													 upperRow))
	{
        // The node may not be strictly within the tile, it can overlap the borders.
        std::unordered_set<CoordValueType> borderPixels;
        node.m_Contour.GenerateBorderPixels(borderPixels, inputLSImageWidth);

        // Loop over the pixels
        for(const auto& pix : borderPixels)
        {
            if(SpatialTools::IsPixelAtTileBorder(pix,
                                                 tile.m_Tx,
                                                 tile.m_Ty,
                                                 lowerCol,
                                                 lowerRow,
                                                 upperCol,
                                                 upperRow,
                                                 nbTilesX,
                                                 nbTilesY,
                                                 inputLSImageWidth))
            {
                listOfBorderNodes.push_back(&node);
                break;
            }            

        } // end for(const auto& pix : borderPixels)

        }

    }; // end of lambdaOp

    graph->ApplyForEachNode(lambdaOp);

    return listOfBorderNodes;
}



template< typename TGraph >
std::unordered_map< typename GraphOperations<TGraph>::NodeType*, uint32_t >
GraphOperations<TGraph>::ExtractStabilityMargin(const GraphPointerType graph,
                        const uint32_t numMaxAdjacencyLayers,
                        const ProcessingTile& tile,
                        const uint32_t nbTilesX,
                        const uint32_t nbTilesY,
                        const uint32_t inputLSImageWidth)
{
    // Retrieve the list of the nodes located at the borders of the tile
    auto borderNodes = GetListOfBorderNodes(graph, tile, nbTilesX, nbTilesY, inputLSImageWidth);
    std::unordered_map< NodeType*, uint32_t > extractedNodes;

    // For each node we visit its adjacency layers using the DFS graph algorithm.
    for(auto& n : borderNodes)
    {
      RecursiveDepthBreadFirstExploration(n, extractedNodes, graph, 0, numMaxAdjacencyLayers);
    }

    return extractedNodes;
}

template <typename TGraph >
void
GraphOperations<TGraph>::RecursiveDepthBreadFirstExploration(NodeType * node,
                                                             std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                                                             const GraphPointerType graph,
                                                             const uint32_t currentNumberOfAdjacencyLayers,
                                                             const uint32_t numMaxAdjacencyLayers)
{
    // First stopping criterion: the max number of adjacency layers is reached.
    if(currentNumberOfAdjacencyLayers > numMaxAdjacencyLayers)
    {
        return;
    }
    else
    {
        // Try to retrieve the depth exploration of this node
        auto isHere = visitedNodes.find(node);

        // The node has been already visited
        if(isHere != visitedNodes.end())
        {

            // If the current depth exploration is lower than the node 's depth exploration
            // then we update its depth with this one and we re-explore this node.
            if(currentNumberOfAdjacencyLayers < isHere->second)
            {
                // Update the depth exploration
                  isHere->second = currentNumberOfAdjacencyLayers;
            
                // Loop over its edges to explore recursively its adjacent nodes.
                for(auto edg : node->m_Edges)
                {
                      RecursiveDepthBreadFirstExploration(graph->GetNodeAt(edg.m_TargetId), 
                                                          visitedNodes, 
                                                          graph, 
                                                          currentNumberOfAdjacencyLayers + 1, 
                                                          numMaxAdjacencyLayers);
                } // end for for(auto edg : node->m_Edges)

              } // end if(currentNumberOfAdjacencyLayers < isHere->second)
            else
            {
                // Second stopping criterion: this node has been fully explored.
                return;
            }
        } // end if
        else
        {
            // This node is explored for the first time.
            visitedNodes[node] = currentNumberOfAdjacencyLayers;

            // Loop over its edges to explore recursively its adjacent nodes.
            for(auto edg : node->m_Edges)
            {

                RecursiveDepthBreadFirstExploration(graph->GetNodeAt(edg.m_TargetId), 
                                                    visitedNodes, 
                                                    graph, 
                                                    currentNumberOfAdjacencyLayers + 1, 
                                                    numMaxAdjacencyLayers);
            } // end for(auto edg : node->m_Edges)

        } // end else

    } // end else

}
/**TODO*/
template <typename TGraph >
template<typename LambdaFunc>
void
GraphOperations<TGraph>
::RecursiveDepthBreadFirstExploration(NodeType * node,
                                      std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                                      const GraphPointerType graph,
                                      const uint32_t currentNumberOfAdjacencyLayers,
                                      LambdaFunc &lambda)
{
    // First stopping criterion according to lambda function
    if(lambda)
    {
        return;
    }
    else
    {
        // Try to retrieve the depth exploration of this node
        auto isHere = visitedNodes.find(node);

        // The node has been already visited
        if(isHere != visitedNodes.end())
        {

            // If the current depth exploration is lower than the node 's depth exploration
            // then we update its depth with this one and we re-explore this node.
            if(currentNumberOfAdjacencyLayers < isHere->second)
            {
                // Update the depth exploration
                isHere->second = currentNumberOfAdjacencyLayers;

                // Loop over its edges to explore recursively its adjacent nodes.
                for(auto edg : node->m_Edges)
                {
                    RecursiveDepthBreadFirstExploration(graph->GetNodeAt(edg.m_TargetId),
                                                        visitedNodes,
                                                        graph,
                                                        currentNumberOfAdjacencyLayers + 1,
                                                        lambda);
                } // end for for(auto edg : node->m_Edges)

            } // end if(currentNumberOfAdjacencyLayers < isHere->second)
            else
            {
                // Second stopping criterion: this node has been fully explored.
                return;
            }
        } // end if
        else
        {
            // This node is explored for the first time.
            visitedNodes[node] = currentNumberOfAdjacencyLayers;

            // Loop over its edges to explore recursively its adjacent nodes.
            for(auto edg : node->m_Edges)
            {

                RecursiveDepthBreadFirstExploration(graph->GetNodeAt(edg.m_TargetId),
                                                    visitedNodes,
                                                    graph,
                                                    currentNumberOfAdjacencyLayers + 1,
                                                    lambda);
            } // end for(auto edg : node->m_Edges)

        } // end else

    } // end else
}


template <typename TGraph >
void
GraphOperations<TGraph>::AggregateGraphs(GraphPointerType graph, 
                     GraphPointerType otherGraph)
{
    // Step 1: increment the id of the nodes in otherGraph by the number of nodes
    // in graph in order to preserve the relative order of the nodes
    // in the adjacency list of the resulting graph.
    const IdType numNodes = graph->GetNumberOfNodes();

    auto lambdaIncrementIdEdge = [&](EdgeType& edge){
        edge.m_TargetId = edge.m_TargetId + numNodes;
    };
    
    auto lambdaOp = [&](NodeType& node)
    {
    	node.m_Id = node.m_Id + numNodes;
    	node.ApplyForEachEdge(lambdaIncrementIdEdge);
    };
    
    otherGraph->ApplyForEachNode(lambdaOp);
    
    // Step 2: We can safely copy the nodes of subgraph into the adjacent list of graph.
    graph->InsertAtEnd(otherGraph);
}


template< typename TGraph >
std::unordered_map<CoordValueType, std::vector< typename GraphOperations<TGraph>::NodeType*> >
GraphOperations<TGraph>::BuildBorderNodesMapForFinalAggregation(const GraphPointerType graph,
                                                                std::unordered_set<uint32_t>& rowBounds,
                                                                std::unordered_set<uint32_t>& colBounds,
                                                                const uint32_t inputLSImageWidth)
{
  // Map : pix coords -> list of nodes.
  std::unordered_map<CoordValueType, std::vector< NodeType*> > borderNodeMap;

    // Loop over the nodes
    for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    {
        // Generate the border pixels
        std::unordered_set< CoordValueType > borderPixels;
        nodeIt->m_Contour.GenerateBorderPixels(borderPixels, inputLSImageWidth);

        // Loop over the border pixels
        for(auto& pix : borderPixels)
        {

            if(SpatialTools::IsPixelAtBoundaries(pix,
                                                 rowBounds,
                                                 colBounds,
                                                 inputLSImageWidth))
            {
                borderNodeMap[pix].push_back(&(*nodeIt));
            }

        } // end for (auto& pix : borderPixels)

    } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)

      return borderNodeMap;
}

template <typename TGraph >
std::unordered_map<CoordValueType, std::vector< typename GraphOperations<TGraph>::NodeType*> > 
GraphOperations<TGraph>::BuildBorderNodesMap(const GraphPointerType graph,
                                             const ProcessingTile& tile,
                                             const uint32_t nbTilesX,
                                             const uint32_t nbTilesY,
                                             const uint32_t inputLSImageWidth)
{
      std::unordered_set< uint32_t > rowBounds;
      std::unordered_set< uint32_t > colBounds;
      std::unordered_map<CoordValueType, std::vector< NodeType*> > borderNodeMap;

    if(tile.m_Ty > 0)
    {
        rowBounds.insert(tile.m_Frame.GetIndex(1));
        rowBounds.insert(tile.m_Frame.GetIndex(1) - 1);
    }

    if(tile.m_Ty < nbTilesY - 1)
    {
        rowBounds.insert(tile.m_Frame.GetIndex(1) + tile.m_Frame.GetSize(1) - 1);
        rowBounds.insert(tile.m_Frame.GetIndex(1) + tile.m_Frame.GetSize(1));
    }

    if(tile.m_Tx > 0)
    {
        colBounds.insert(tile.m_Frame.GetIndex(0));
        colBounds.insert(tile.m_Frame.GetIndex(0) - 1);
    }

    if(tile.m_Tx < nbTilesX - 1)
    {
        colBounds.insert(tile.m_Frame.GetIndex(0) + tile.m_Frame.GetSize(0) - 1);
        colBounds.insert(tile.m_Frame.GetIndex(0) + tile.m_Frame.GetSize(0));
    }

    for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
    {
        // Generate the border pixels
        std::unordered_set< CoordValueType > borderPixels;
        nodeIt->m_Contour.GenerateBorderPixels(borderPixels, inputLSImageWidth);

        // Loop over the border pixels
        for(auto& pix : borderPixels)
        {

            if(SpatialTools::IsPixelAtBoundaries(pix,
                                                 rowBounds,
                                                 colBounds,
                                                 inputLSImageWidth))
            {
                borderNodeMap[pix].push_back(&(*nodeIt));
            }

        } // end for (auto& pix : borderPixels)

    } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)

    return borderNodeMap;
}

template <typename TGraph >
void
GraphOperations<TGraph>
::RemoveDuplicatedNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                        GraphPointerType graph,
                        const uint32_t inputLSImageWidth)
{

    std::unordered_set< CoordValueType > pixelsToRemove;

    // Loop over the pixels on the common tile borders.
    for(auto& kv : borderNodeMap)
    {
        // If this pixel belongs to a node already processed or if there are no
        // duplicated nodes then we do nothing.
        if(pixelsToRemove.find(kv.first) == pixelsToRemove.end() && kv.second.size() > 1)
        {
            // First step: sort duplicated nodes wrt to their id in order to preserve
              // the stability of this algorithm and the relative order of the nodes
              // between the tiles.
              std::sort(kv.second.begin(), kv.second.end(), [](const NodeType * a, 
                                                              const NodeType * b)->bool{
                return a->m_Id < b->m_Id;    
            });

              // Step 1: determine the node which will stay (the one with the lower id)
              auto remainingNode = kv.second.front();
              //uint64_t refStartingCoords = remainingNode->GetFirstPixelCoords();

              // Step 2: Loop over the duplicated nodes
              auto nodePtrIt = kv.second.begin();
              nodePtrIt++;

            for(; nodePtrIt != kv.second.end(); nodePtrIt++)
            {
                // Loop over their edges
                for(auto& edg : (*nodePtrIt)->m_Edges)
                {
                    // Retrieve the adjacent node to the duplicated node to remove.
                      auto adjToDNode = graph->GetNodeAt(edg.m_TargetId);

                      // Find the edge targeting the duplicated node
                      auto edgItAdjToDNode = adjToDNode->FindEdge((*nodePtrIt)->m_Id);

                      // Keep in memory the boundary value if an edge has to be created
                      // for the reference node.
                      uint32_t boundary = edgItAdjToDNode->m_Boundary;

                      // Determine if the adjacent node of the duplicated node is an adjacent
                      // node of the reference node.
                      auto isAdjToRefNode = std::find_if(remainingNode->m_Edges.begin(), remainingNode->m_Edges.end(),
                                           [&](const EdgeType& edgIn)->bool{
                                               return (graph->GetNodeAt(edgIn.m_TargetId))->GetFirstPixelCoords() == adjToDNode->GetFirstPixelCoords();
                    });

                    // Determine according to the previous operation if edges has to be created.
                    bool createEdge = (isAdjToRefNode == remainingNode->m_Edges.end() ) ? true : false;

                    // Can safely remove the edge between the duplicated node and its adjacent node.
                    adjToDNode->m_Edges.erase(edgItAdjToDNode);

                    if(createEdge)
                    {
                        // Create an edge: refNode -> adjToDNode
                        auto refToAdj = remainingNode->AddEdge();
                        refToAdj->m_TargetId = adjToDNode->m_Id;
                        refToAdj->m_Boundary = boundary;

                        // Create an edge: adjToDNode -> refNode
                        auto adjToRef = adjToDNode->AddEdge();
                        adjToRef->m_TargetId = remainingNode->m_Id;
                        adjToRef->m_Boundary = boundary;

                    } // end if(isAdjToRefNode == adjToDNode->m_Edges.end())        

                } // end for(auto& edg : (*nodePtrIt)->m_Edges)

                  (*nodePtrIt)->m_HasToBeRemoved = true;

            } // end for(; nodePtrIt != kv.second.end(); nodePtrIt++)

            // Remove the pixels of those duplicated nodes (they must all have the same set of pixels)
            std::unordered_set< CoordValueType > borderPixels;
            remainingNode->m_Contour.GenerateBorderPixels(borderPixels, inputLSImageWidth);
            for(const auto& pix : borderPixels)
            {

                pixelsToRemove.insert(pix);
    
            } // end for(const auto& pix : borderPixels)

        } // end if(kv.second.size() > 1)

    } // end for(auto& kv : borderNodeMap)

    // Remove the pixels of the duplicated node from the
    // map.
    for(const auto& pix : pixelsToRemove)
    {
        // TODO : cppcheck tells us this check is useless, remove?
        if(borderNodeMap.find(pix) != borderNodeMap.end())
        {
              borderNodeMap.erase(pix);
        }
    }
}

template <typename TGraph >
void 
GraphOperations<TGraph>::DetectNewAdjacentNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                                                GraphPointerType graph,
                                                const uint32_t inputLSImageWidth,
                                                const uint32_t inputLSImageHeight)
{
    // Loop over the pixels that are within just one node.
    for(auto& kv : borderNodeMap)
    {
          // Determine the 4-neighborhood of the pixel
          auto neighbors = SpatialTools::FourConnectivity(kv.first, inputLSImageWidth, inputLSImageHeight);

          // Retrieve the node that contains this pixel
          auto currBorderNode = kv.second.front();


          // Loop over the neigbhoring pixels.
          for(unsigned short n = 0; n < 4; n++)
        {
            // Is it a pixel valid ?
              if(neighbors[n] > -1)
            {

                  // Look if the neighbor pixel is in the map
                  auto isNeighInMap = borderNodeMap.find(neighbors[n]);

                  if(isNeighInMap != borderNodeMap.end())
                {

                      // Retrieve the node that contains the neighboring pixel
                      auto neighBorderNode = isNeighInMap->second.front();

                      // If neighBorderPixel is different than currBorderNode
                      // then it may be a new adjacent node.
                      if(neighBorderNode != currBorderNode)
                    {

                          // Determine if neighBorderNode is already an adjacent
                          // node of currBorderNode
                          auto currToNeigh = currBorderNode->FindEdge(neighBorderNode->m_Id);


                          if(currToNeigh == currBorderNode->m_Edges.end())
                        {
                              // neighBorderNode is a new adjacent node,
                              // the first step consists of determining the
                              // boundary length between currBorderNode and
                              // neighBorderNode.
                              uint32_t boundary = 0;

                              // Generate the border pixels of currBorderNode
                              std::unordered_set< CoordValueType > currBorderPixels;
                              currBorderNode->m_Contour.GenerateBorderPixels(currBorderPixels, inputLSImageWidth);

                              // Generate the border pixels of neighBorderNode
                              std::unordered_set< CoordValueType > neighBorderPixels;
                              neighBorderNode->m_Contour.GenerateBorderPixels(neighBorderPixels, inputLSImageWidth);

                              // Loop over the border pixels of currBorderNode
                              for(const auto& pix : currBorderPixels)
                            {

                                  auto neighPixels = SpatialTools::FourConnectivity(pix, inputLSImageWidth, inputLSImageHeight);

                                  for(unsigned short nn = 0; nn < 4; nn++)
                                {

                                      if(neighPixels[nn] > -1 && 
                                         neighBorderPixels.find(neighPixels[nn]) != neighBorderPixels.end())
                                    {
                                          boundary++;
                                    }
                                } // end for(unsigned short nn = 0; nn < 4; nn++)

                            } // end for(const auto& pix : currBorderPixels)

                              // The boundary must not be null
                              assert(boundary>=1);

                              // Step 2: Add edges

                              // Edge : currBorderNode -> neighBorderNode
                              auto currToNeighIn = currBorderNode->AddEdge();
                              currToNeighIn->m_TargetId = neighBorderNode->m_Id;
                              currToNeighIn->m_Boundary = boundary;

                              // Edge : neighBorderNode -> currBorderNode
                              auto neighToCurr = neighBorderNode->AddEdge();
                              neighToCurr->m_TargetId = currBorderNode->m_Id;
                              neighToCurr->m_Boundary = boundary;

                        } // end if(currToNeigh == currBorderNode->m_Edges.end())

                    } // end if(neighBorderNode != currBorderNode)

                } // end if(isNeighInMap != borderNodeMap.end())

            } // end if(neighbors[n] > -1)

        } // end for(unsigned short n = 0; n < 4; n++)

    } // end for(auto& kv : borderNodeMap)

  // Can safely remove the duplicated nodes.
  graph->RemoveNodes();
}
template <typename TGraph >
void GraphOperations<TGraph>
::DisplayNode(const GraphPointerType graph, unsigned int coordNode)
 {
	std::cout << "Looking for node " << coordNode << std::endl;
	for(int i = 0; i < graph->GetNumberOfNodes(); ++i)
	{
		auto node = graph->GetNodeAt(i);
		if(node->GetFirstPixelCoords() == coordNode)
		{
			for(int k = 0; k < node->m_Edges.size(); k++)
			{
				auto adjNode =  graph->GetNodeAt(node->m_Edges[k].m_TargetId);
				std::cout << "PROC " << MPIConfig::Instance()->GetMyRank() << " :: Adjacent Marge stabilité = " << adjNode->GetFirstPixelCoords() << std::endl;
			}
		}
	}

 }

} // end of namespace obia
} // end of namespace otb

#endif

