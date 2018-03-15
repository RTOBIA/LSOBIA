#ifndef otbObiaGraphFileReader_txx
#define otbObiaGraphFileReader_txx
#include "otbSystem.h"
#include <itksys/SystemTools.hxx>
#include "otbObiaGraphFileReader.h"
#include "otbObiaGraphOperations.h"
#include "itkIntTypes.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <unordered_map>

namespace otb
{
namespace obia
{

template< typename TOutputGraph >
GraphFileReader< TOutputGraph >
::GraphFileReader()
 :m_FileName("")
{
    // Create the output. We use static_cast<> here because we know the default
    // output must be of type TOutputGraph
    typename TOutputGraph::Pointer output =
    static_cast< TOutputGraph * >( this->MakeOutput(0).GetPointer() );
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput( 0, output.GetPointer() );

    // Set the default behavior of an image source to NOT release its
    // output bulk data prior to GenerateData() in case that bulk data
    // can be reused (an thus avoid a costly deallocate/allocate cycle).
    this->ReleaseDataBeforeUpdateFlagOff();
}

template<typename TOutputGraph>
void
GraphFileReader< TOutputGraph >::
PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "m_FileName: " << this->m_FileName << "\n";
}

template<typename TOutputGraph>
void
GraphFileReader< TOutputGraph >::
UpdateOutputInformation()
{
    // Mise a jour du graphe en fonction de son header
    typename TOutputGraph::Pointer graph = this->GetOutput();
    GraphImageInfo info = GraphOperations<TOutputGraph>::ReadGraphHeader(m_FileName);
    graph->SetImageWidth(info.m_width);
    graph->SetImageHeight(info.m_height);
    graph->SetNumberOfSpectralBands(info.m_nbBands);
    graph->SetProjectionRef(info.m_projectionRef);

    itk::ModifiedTimeType t1;
    t1 = this->GetMTime();
    graph->SetPipelineMTime(t1);

}


template< typename TOutputGraph >
void GraphFileReader< TOutputGraph >
::GenerateData()
{
    typename TOutputGraph::Pointer graph = this->GetOutput();

    // Open the file stream
    std::ifstream inFile(m_FileName, std::ios::in | std::ios::binary);
    if(inFile.good())
    {
        // Retrieve the number of bytes to read.
        uint64_t numberOfBytes;
        inFile.read(reinterpret_cast<char*>(&numberOfBytes), UInt64Size);

        // Read the file
        std::vector<char> serializedGraph(numberOfBytes);
        inFile.read(&serializedGraph[0], numberOfBytes * CharSize);
        inFile.close();

        // Local variable: offset in the bit stream.
        uint64_t position = 0;

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
            if(startPixIdMap.find(newNode->GetFirstPixelCoords()) == startPixIdMap.end())
            {
                startPixIdMap[newNode->GetFirstPixelCoords()] = id;
            }
            else
            {
                std::cerr << "Error: there are least two nodes sharing the same starting coordinates." << std::endl;
                exit(EXIT_FAILURE);
            }

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
		from_stream(serializedGraph,newEdge->m_TargetId, position);

                // The boundary
		from_stream(serializedGraph,newEdge->m_Boundary,position);

                // Serialize the specific attributes of the edges
                newEdge->m_Attributes.DeSerialize(serializedGraph, position);
            }
        }

        // Update the target ids of the edges with the real ids of the nodes
        for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)
        {
            for(auto& edg : nodeIt->m_Edges)
            {
                auto isInStartPixIdMap = startPixIdMap.find(edg.m_TargetId);

                if(isInStartPixIdMap == startPixIdMap.end())
                {
                    std::cerr << "Error: the starting coordinates of the edge is not in the map." << std::endl;
                    exit(EXIT_FAILURE);
                }
                else
                {
                    edg.m_TargetId = isInStartPixIdMap->second;
                }
            }
        }
    }
    else
    {
        std::cerr << "Cannot open input file to load the graph: " << m_FileName << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename TOutputGraph >
void GraphFileReader< TOutputGraph >
::TestFileExistanceAndReadability()
{
    // Test if the file exists.
    if (!itksys::SystemTools::FileExists(this->m_FileName.c_str()))
    {
        otb::obia::GraphFileReaderException e(__FILE__, __LINE__);
        std::ostringstream msg;
        msg << "The file doesn't exist. "
            << std::endl << "Filename = " << this->m_FileName
            << std::endl;
        e.SetDescription(msg.str().c_str());
        throw e;
    }

    // Test if the file can be open for reading access.
    //Only if m_FileName specify a filename (not a dirname)
    if (itksys::SystemTools::FileExists(this->m_FileName.c_str(), true) == true)
    {
        std::ifstream readTester;
        readTester.open(this->m_FileName.c_str());
        if (readTester.fail())
        {
            readTester.close();
            std::ostringstream msg;
            msg << "The file couldn't be opened for reading. "
                << std::endl << "Filename: " << this->m_FileName
                << std::endl;
            otb::obia::GraphFileReaderException e(__FILE__, __LINE__, msg.str().c_str(), ITK_LOCATION);
            throw e;
        }
        readTester.close();
    }
}

} // end of namespace obia
} // end of namespace otb

#endif
