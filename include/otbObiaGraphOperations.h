/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef otbObiaGraphOperations_h
#define otbObiaGraphOperations_h
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include "otbObiaGraph.h"
#include <unordered_set>

/**
\file otbObiaGraphOperations.h
\brief This file define the class used for several operations in graph like extracting staiblity margins, serializing ...
*/
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
        double m_originX;
        double m_originY;

        GraphImageInfo():
            m_width(0),
            m_height(0),
            m_nbBands(0),
            m_projectionRef(""),
			m_originX(0),
			m_originY(0)
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

      /**\brief  This method removes the nodes which are fully within the margin area.
       *  Those nodesmight be instable. Therefore, in order to avoid propagating errors in the obia
       *  chains, those nodes have to be removed from the adjacent graph. The referential
       *  used to detect the unstable segment is the large scale input image referential.
       *  It is assumed that the starting coordinates and the bounding boxes of the nodes
       *  are in this referential. Otherwise, the resulting graph is unstable.\n
       *  \param graph: adjacent graph
       *  \param imageWidth: input large scale image width*/
      static void RemoveUnstableNodes(GraphPointerType graph,
                      const ProcessingTile& tile,
                      const uint32_t imageWidth);

      /**\brief This method serializes an extracted stability margin into a sequence of bits that
       can be either sent via MPI requests or written in a binary file.
       \param: Stability margin
       \param: Graph used to check if a node contains a target which is not in the stability margins*/
      static std::vector<char> SerializeStabilityMargin(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                            const GraphPointerType graph);

      /**\brief This method serializes an adjacent graph into a sequence of bits that
       *  can be either sent via MPI requests or written in a binary file.
       * \param: Graph to serialize*/
      static std::vector<char> SerializeGraph(const GraphPointerType graph);

       /**\brief This methods builds a graph from a bit stream. It is generic if all the needed
        * serialization methods of the specific attributes are provided by the user.
        * \param: Deserialized graph
        * \return: Pointer to reconstruced graph*/
      static GraphPointerType DeSerializeGraph(const std::vector<char>& serializedGraph);

      /**\brief This methods serialize and write the margin graph margin to disk
       * \param : Stability margin
       * \param : Graph pointer used to serialize stability margins
       * \param : File path for output file*/
      static void WriteMarginGraphToDisk(std::unordered_map< NodeType*, uint32_t >& stabilityMargin,
                     const GraphPointerType graph,
                     const std::string& outputPath);

      /**\brief This methods writes the serialized stability margin to disk.
       * \param: Stability margin
       * \param: File path for output file*/
      static void WriteSerializedMarginToDisk(const std::vector<char>& serializedStabilityMargin,
                          const std::string& outputPath);

      /**\brief This methods simply writes the bit stream to a file stored on the disk.
       * \param: Graph pointer
       * \param: Output file path */
      static void WriteGraphToDisk(const GraphPointerType graph,
                   const std::string& outputPath);

      /**\brief This methods loads a graph from its storage file on the disk.
       * \param: Input file path
       * \return: Pointer graph*/
      static GraphPointerType ReadGraphFromDisk(const std::string & inputPath);

      /**\brief This methods loads a serialized margin from its storage file on the disk.
      * \param: Input file path
      * \return: Vector of char corresponding to the serialized margin*/
      static std::vector<char> ReadSerializedMarginFromDisk(const std::string& inputPath);

      /**\brief Give the tile information (dimension and location in the input image), this method
       * returns a list of nodes located at the borders of the tiles.
       * \param: Graph which border nodes will be extracted
       * \param: Current tile processed
       * \param: Number of tiles X
       * \param: Number of tiles Y
       * \param: Width of the image (used to check if a node is border)
       * \return: A list of border nodes*/
      //
      static std::vector< NodeType* > GetListOfBorderNodes(const GraphPointerType graph,
                                                           const ProcessingTile& tile,
                                                           const uint32_t nbTilesX,
                                                           const uint32_t nbTilesY,
                                                           const uint32_t inputLSImageWidth);

      /**\brief This methods extracts a subgraph corresponding to the stability margin for graphs of
       * adjacent tiles. The stability margin value corresponds to the maximum number of adjacency layers.
       * \param: Graph used to extract stability margin
       * \param: Number of adjacency layers
       * \param: Processed tile
       * \param: Number of tiles X
       * \param: Number of tiles Y
       * \param: Width of the image
       * \return: A list of nodes corresponding to the stability margin*/
      static std::unordered_map< NodeType*, uint32_t > ExtractStabilityMargin(const GraphPointerType graph,
                                          const uint32_t numMaxAdjacencyLayers,
                                          const ProcessingTile& tile,
                                          const uint32_t nbTilesX,
                                          const uint32_t nbTilesY,
                                          const uint32_t inputLSImageWidth);

      /**\brief This methods consists of a depth-bread first exploration of the edges of node given
       * the map of the visited nodes so far with their exploration depth, the current exploration
       * depth and the maximum exploration depth.
       * \param: Node to explore
       * \param: List of visited nodes
       * \param: Graph to get nodes
       * \param: Current number of adjacency layers
       * \param: Max adjacency layers number*/

      static void RecursiveDepthBreadFirstExploration(NodeType * node,
                              std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                              const GraphPointerType graph,
                              const uint32_t currentNumberOfAdjacencyLayers,
                              const uint32_t numMaxAdjacencyLayers);

      /**\brief This methods consists of a depth-bread first exploration of the edges of node given
       * the map of the visited nodes so far with their exploration depth, the current exploration
       * depth and the maximum exploration depth.
       * \param: Node to explore
       * \param: List of visited nodes
       * \param: Graph to get nodes
       * \param: Current number of adjacency layers
       * \param: Lambda function*/
      template<typename LambdaFunc>
      static void RecursiveDepthBreadFirstExploration(NodeType * node,
                                                      std::unordered_map< NodeType*, uint32_t >& visitedNodes,
                                                      const GraphPointerType graph,
                                                      const uint32_t currentNumberOfAdjacencyLayers,
                                                      LambdaFunc &lambda);

      /**\brief This methods aggregates two graphs while preserving the property of the node ids.
       * \param: Graph
       * \param: Graph to insert*/
      static void AggregateGraphs(GraphPointerType graph, 
                  GraphPointerType otherGraph);


      /**\brief This method builds a map that assigns for each pixel located on the common borders
       * of the tile with the adjacent tiles a node
       * \param: Graph
       * \param: Processed tile
       * \param: Number of tiles X
       * \param: Number of tiles Y
       * \param: image width
       * \return map of borders coordinates and vector of node associated to coordinate
       * */
      static std::unordered_map<CoordValueType, std::vector<NodeType*> > BuildBorderNodesMap(const GraphPointerType graph,
                                                 const ProcessingTile& tile,
                                                 const uint32_t nbTilesX,
                                                 const uint32_t nbTilesY,
                                                 const uint32_t inputLSImageWidth);

      /**\brief This methods does the same thing than the method BuildBorderNodesMap. The only difference
       * is that boundaries are given as parameters. This may be be factorized for a cleaner code.
       * \param: Graph
       * \param: row boundaries
       * \parma: col boundaries
       * \param: image width
       * \return map of borders coordinates and vector of node associated to coordinate
       * */
      static std::unordered_map<CoordValueType, std::vector<NodeType*> > BuildBorderNodesMapForFinalAggregation(const GraphPointerType graph,
                                                        std::unordered_set<uint32_t>& rowBounds,
                                                        std::unordered_set<uint32_t>& colBounds,
                                                        const uint32_t inputLSImageWidth);

      /**\brief This methods removes the duplicated nodes from the graph.
       * \param: Map of border node (only border node can be duplicated)
       * \param: Graph containing nodes to remove
       * \param: Image width*/
      static void RemoveDuplicatedNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                    GraphPointerType graph,
                    const uint32_t inputLSImageWidth);


      /**\brief This method detects the new pair of adjacent nodes on the common tile borders.
       * \param: Border nodes map
       * \param: Graph
       * \param: Image width
       * \param: Image height*/
      static void DetectNewAdjacentNodes(std::unordered_map<CoordValueType, std::vector< NodeType*> >& borderNodeMap,
                     GraphPointerType graph,
                     const uint32_t inputLSImageWidth,
                     const uint32_t inputLSImageHeight);

      /**\brief This method write a graph header containing meta data used for process like image width
       * \param: File path
       * \param: Graph containing metadata
       */
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
            content << graph->GetOriginX() <<"\n";
            content << graph->GetOriginY() <<"\n";
            content << graph->GetProjectionRef()<<"\n";
            std::ofstream header;
            header.open(headerPath);
            header << content.str();
            header.close();
      }



      /**\brief This method read a graph header containing meta data used for process like image width
      * \param: File path
      * \return: Graph information
      */
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

            //Read height
            std::getline (header,line);
            info.m_height = std::stod(line);

            //Read width
            std::getline (header,line);
            info.m_width = std::stod(line);

            //Read nb bands
            std::getline (header,line);
            info.m_nbBands = std::stod(line);

            //Read originX
            std::getline (header,line);
            info.m_originX = std::stod(line);

            //Read originY
			std::getline (header,line);
		    info.m_originY = std::stod(line);

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

    /**\brief This method display node and its information (used for debug)
     * \param: Graph
     * \param: coordinate of node to display
     */
    static void DisplayNode(const GraphPointerType graph, unsigned int coordNode);

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphOperations.txx"
#endif
