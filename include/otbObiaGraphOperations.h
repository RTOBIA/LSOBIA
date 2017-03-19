#ifndef otbObiaGraphOperations_h
#define otbObiaGraphOperations_h
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include "otbObiaGraph.h"

namespace otb
{
namespace obia
{

    /* Helper class */
    class GraphImageInfo{
    public:
        int m_width;
        int m_height;
        int m_nbBands;
        std::string m_projectionRef;

        GraphImageInfo():
            m_width(0),
            m_height(0),
            m_nbBands(0),
            m_projectionRef("")
        {}
    };
    /*
     *    LSGraphOperations
     *
     *    This class is like a toolbox
     *    containing graph operations for
     *    large scale processes (sequential tiling and
     *    distributed execution)
     **/
    template< typename TGraph >
    class GraphOperations
    {

      public:

      /* Some convenient alias */
      using GraphType = TGraph;
      using GraphPointerType = typename GraphType::Pointer;
      using NodeType = typename GraphType::NodeType;
      using EdgeType = typename GraphType::EdgeType;

      // This method removes the nodes which are fully within the margin area. Those nodes
      // might be instable. Therefore, in order to avoid propagating errors in the obia
      // chains, those nodes have to be removed from the adjacent graph. The referential
      // used to detect the unstable segment is the large scale input image referential.
      // It is assumed that the starting coordinates and the bounding boxes of the nodes
      // are in this referential. Otherwise, the resulting graph is unstable.
      //
      // @param graph: adjacent graph
      //    @param imageWidth: input large scale image width

      static void RemoveUnstableNodes(GraphPointerType graph,
                      const ProcessingTile& tile,
                      const uint32_t imageWidth);

      // This method serializes an extracted stability margin into a sequence of bits that
      // can be either sent via MPI requests or written in a binary file.
      static std::vector<char> SerializeStabilityMargin(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                            const GraphPointerType graph);

      // This method serializes an adjacent graph into a sequence of bits that
      // can be either sent via MPI requests or written in a binary file.
      static std::vector<char> SerializeGraph(const GraphPointerType graph);

      // This methods builds a graph from a bit stream. It is generic if all the needed
      // serialization methods of the specific attributes are provided by the user.
      static GraphPointerType DeSerializeGraph(const std::vector<char>& serializedGraph);

      // This methods serialize and write the serialized stability margin to disk.
      static void WriteMarginGraphToDisk(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                     const GraphPointerType graph,
                     const std::string& outputPath);

      // This methods writes the serialized stability margin to disk.
      static void WriteSerializedMarginToDisk(const std::vector<char>& serializedStabilityMargin,
                          const std::string& outputPath);

      // This methods simply writes the bit stream to a file stored on the disk.
      static void WriteGraphToDisk(const GraphPointerType graph,
                   const std::string& outputPath);

      // This methods loads a graph from its storage file on the disk.
      static GraphPointerType ReadGraphFromDisk(const std::string & inputPath);

      static std::vector<char> ReadSerializedMarginFromDisk(const std::string& inputPath);

      // Give the tile information (dimension and location in the input image), this method
      // returns a list of nodes located at the borders of the tiles.
      static std::vector< NodeType* > GetListOfBorderNodes(const GraphPointerType graph,
                                                           const ProcessingTile& tile,
                                                           const uint32_t nbTilesX,
                                                           const uint32_t nbTilesY,
                                                           const uint32_t inputLSImageWidth);

      // This methods extracts a subgraph corresponding to the stability margin for graphs of
      // adjacent tiles. The stability margin value corresponds to the maximum number of adjacency
      // layers.
      static std::unordered_map< NodeType*, uint32_t > ExtractStabilityMargin(const GraphPointerType graph,
                                          const uint32_t numMaxAdjacencyLayers,
                                          const ProcessingTile& tile,
                                          const uint32_t nbTilesX,
                                          const uint32_t nbTilesY,
                                          const uint32_t inputLSImageWidth);

      // This methods consists of a depth-bread first exploration of the edges of node given
      // the map of the visited nodes so far with their exploration depth, the current exploration
      // depth and the maximum exploration depth.
      static void RecursiveDepthBreadFirstExploration(NodeType * node,
                              std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                              const GraphPointerType graph,
                              const uint32_t currentNumberOfAdjacencyLayers,
                              const uint32_t numMaxAdjacencyLayers);

      // This methods consists of a depth-bread first exploration of the edges of node given
      // the map of the visited nodes so far with their exploration depth, the current exploration
      // depth and the maximum exploration depth.
      template<typename LambdaFunc>
      static void RecursiveDepthBreadFirstExploration(NodeType * node,
                                                                          std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                                                      const GraphPointerType graph,
                                                      const uint32_t currentNumberOfAdjacencyLayers,
                                                      LambdaFunc &lambda);

      // This methods aggregates two graphs while preserving the property of the node ids.
      static void AggregateGraphs(GraphPointerType graph, 
                  GraphPointerType otherGraph);


      // This method builds a map that assigns for each pixel located on the common borders
      // of the tile with the adjacent tiles a node
      static std::unordered_map<CoordValueType, std::vector<NodeType*> > BuildBorderNodesMap(const GraphPointerType graph,
                                                 const ProcessingTile& tile,
                                                 const uint32_t nbTilesX,
                                                 const uint32_t nbTilesY,
                                                 const uint32_t inputLSImageWidth);

      // This methods does the same thing than the method BuildBorderNodesMap. The only difference 
      // is that boundaries are given as parameters. This may be be factorized for a cleaner code.
      static std::unordered_map<CoordValueType, std::vector<NodeType*> > BuildBorderNodesMapForFinalAggregation(const GraphPointerType graph,
                                                        std::unordered_set<uint32_t>& rowBounds,
                                                        std::unordered_set<uint32_t>& colBounds,
                                                        const uint32_t inputLSImageWidth);

      // This methods removes the duplicated nodes from the graph.
      static void RemoveDuplicatedNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                    GraphPointerType graph,
                    const uint32_t inputLSImageWidth);


      // This method detects the new pair of adjacent nodes on the common tile borders.
      static void DetectNewAdjacentNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                     GraphPointerType graph,
                     const uint32_t inputLSImageWidth,
                     const uint32_t inputLSImageHeight);

      static void WriteGraphHeader(const std::string & outputPath, GraphPointerType graph)
      {
            // Write header file (contains metadata)
            std::string headerPath(outputPath);
            headerPath+=".hdr";

            std::ostringstream content;
            content << outputPath<<"\n";
            content << graph->GetImageHeight()<<"\n";
            content << graph->GetImageWidth()<<"\n";
            content << graph->GetNumberOfSpectralBands()<<"\n";
            content << graph->GetProjectionRef()<<"\n";

            std::ofstream header;
            header.open(headerPath);
            header << content.str();
            header.close();
      }


    static GraphImageInfo ReadGraphHeader(const std::string & inputPath)
    {
        std::string inputPathHeader = inputPath + ".hdr";
        GraphImageInfo info;
        std::string line;
        std::ifstream header (inputPathHeader);
        if (header.is_open())
        {
            // Skip first line
            std::getline (header,line);
            std::getline (header,line);
            info.m_height = std::stod(line);
            std::getline (header,line);
            info.m_width = std::stod(line);
            std::getline (header,line);
            info.m_nbBands = std::stod(line);
            std::string projectionRef;
            while (std::getline (header,line))
            {
                projectionRef+=line+"\n";
            }
            info.m_projectionRef = projectionRef;
            header.close();
            return info;
        }else
        {
            return GraphImageInfo();
        }
    }

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphOperations.txx"
#endif
