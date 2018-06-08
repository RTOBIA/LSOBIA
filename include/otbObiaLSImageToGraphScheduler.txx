#ifndef otbObiaLSImageToGraphScheduler_txx
#define otbObiaLSImageToGraphScheduler_txx
#include "otbObiaLSImageToGraphScheduler.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"

#include "otbObiaStreamingGraphToImageFilter.h"

namespace otb
{
namespace obia
{

template< typename TInputImage, typename TOutputGraph>
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::LSImageToGraphScheduler() :
 m_OutputDir(""),
m_LabelImageName("labelOutput"),
 m_WriteLabelImage(false),
 m_SplittedGraph(false),
 m_WriteGraph(false),
 m_PaddingValue(0),
 m_Graph(nullptr),
 m_MaxTileSizeX(0),
 m_MaxTileSizeY(0),
 m_ImageHeight(0),
 m_ImageWidth(0),
 m_NumberOfTilesX(0),
 m_NumberOfTilesY(0),
 m_MaxNumberOfTilesPerProcessor(0),
 m_OriginX(0),
 m_OriginY(0),
 m_NumberOfSpectralBands(0),
 m_AvailableMemory(0),
 m_ProcessNoData(false),
 m_NoDataValue(0)
 {
 }

template< typename TInputImage, typename TOutputGraph>
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::~LSImageToGraphScheduler()
{

}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::Update()
 {
    this->PreProcessing();
    this->GenerateData();
    this->CreateOutput();
 }

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetFileName(const std::string& filename)
{
    m_FileName = filename;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetLabelImageName(const std::string& labelImage)
{
    m_LabelImageName = labelImage;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetOutputDir(const std::string & path)
{
    m_OutputDir = path;
    if( !m_OutputDir.empty() && m_OutputDir[m_OutputDir.size() - 1] != '/')
    {
        m_OutputDir.append("/");
    }
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetMaxTileSizeX(const uint32_t MaxTileSizeX)
{
    m_MaxTileSizeX = MaxTileSizeX;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetMaxTileSizeY(const uint32_t MaxTileSizeY)
{
    m_MaxTileSizeY = MaxTileSizeY;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetAvailableMemory(const uint64_t mem)
{
    m_AvailableMemory = mem * 1024 * 1024;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetNoDataValue(const float n)
{
    m_NoDataValue = n;
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::SetTemporaryDirectory(const std::string& tmpDir)
{
    m_TemporaryDirectory = tmpDir;

    if(!m_TemporaryDirectory.empty())
    {
        if(m_TemporaryDirectory[m_TemporaryDirectory.size()-1] != '/')
        {
            m_TemporaryDirectory.append("/");
        }
    }
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::PreProcessing()
{
    // Retrieve the input image metadata without allocating it
    auto imgReader = InputImageReaderType::New();
    imgReader->SetFileName(m_FileName);
    imgReader->UpdateOutputInformation();
    m_ImageWidth  = imgReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
    m_ImageHeight = imgReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
    m_NumberOfSpectralBands = imgReader->GetOutput()->GetNumberOfComponentsPerPixel();
    m_ProjectionRef = imgReader->GetOutput()->GetProjectionRef();
    m_OriginX = imgReader->GetOutput()->GetOrigin()[0];
    m_OriginY = imgReader->GetOutput()->GetOrigin()[1];

    // Compute the padding value
    ComputePaddingValue();

    // Re-compute the tile sizes to optimize the load balancing
    ComputeFinalTileSize();

    // Compute its own tiles
    ComputeTileFrames();

    // Determine the maximum number of tiles per processor
    m_MaxNumberOfTilesPerProcessor = 0;
    for(auto& kv : m_TilesPerProcessor)
    {
        if(kv.second.size() > m_MaxNumberOfTilesPerProcessor)
        {
            m_MaxNumberOfTilesPerProcessor = kv.second.size();
        }
    }

    // Create an instance of the graph
    m_Graph = OutputGraphType::New();

}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::ComputeTileFrames()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    uint32_t tid = 0;
    long int rectifiedIndex;
    long int rectifiedPadding;

    for(uint32_t ty = 0; ty < m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < m_NumberOfTilesX; tx++)
        {

            m_TilesPerProcessor[mpiTools->GetProcessorRankFromTileId(tid)].insert(tid);

            if(mpiTools->IsMyTurn(tid))
            {

                ProcessingTile tile;

                tile.m_Tx = tx;
                tile.m_Ty = ty;
                tile.m_Frame.SetIndex(0, tx * m_MaxTileSizeX);
                tile.m_Frame.SetIndex(1, ty * m_MaxTileSizeY);
                tile.m_Frame.SetSize(0, std::min((long int)m_MaxTileSizeX, (long int)m_ImageWidth - tile.m_Frame.GetIndex(0)));
                tile.m_Frame.SetSize(1, std::min((long int)m_MaxTileSizeY, (long int)m_ImageHeight - tile.m_Frame.GetIndex(1)));
                tile.m_MarginValues.fill(0);

                if(ty > 0)
                {
                    rectifiedIndex = std::max((long int)0, (long int)(tile.m_Frame.GetIndex(1) - m_PaddingValue));
                    rectifiedPadding = tile.m_Frame.GetIndex(1) - rectifiedIndex;

                    tile.m_Frame.SetIndex(1, rectifiedIndex);
                    tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) + rectifiedPadding);
                    tile.m_MarginValues[0] = rectifiedPadding;
                }

                if(tx < (m_NumberOfTilesX - 1) )
                {
                    rectifiedPadding = m_PaddingValue;
                    if(tile.m_Frame.GetIndex(0) + tile.m_Frame.GetSize(0) + m_PaddingValue > m_ImageWidth)
                    {
                        rectifiedPadding = m_ImageWidth - (tile.m_Frame.GetIndex(0) + tile.m_Frame.GetSize(0));
                    }

                    tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) + rectifiedPadding);
                    tile.m_MarginValues[1] = rectifiedPadding;
                }

                if(ty < m_NumberOfTilesY - 1)
                {
                    rectifiedPadding = m_PaddingValue;
                    if(tile.m_Frame.GetIndex(1) + tile.m_Frame.GetSize(1) + m_PaddingValue > m_ImageHeight)
                    {
                        rectifiedPadding = m_ImageHeight - (tile.m_Frame.GetIndex(1) + tile.m_Frame.GetSize(1));
                    }

                    tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) + rectifiedPadding);
                    tile.m_MarginValues[2] = rectifiedPadding;
                }

                if(tx > 0)
                {

                    rectifiedIndex = std::max((long int)0, (long int)(tile.m_Frame.GetIndex(0) - m_PaddingValue));
                    rectifiedPadding = tile.m_Frame.GetIndex(0) - rectifiedIndex;

                    tile.m_Frame.SetIndex(0, rectifiedIndex);
                    tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) + rectifiedPadding);
                    tile.m_MarginValues[3] = rectifiedPadding;
                }

                m_TileMap[tid] = tile;

            } // end if (mpiTools->IsMyTurn(tid))

                tid++;

        } // end for(uint32_t tx = 0; tx < nbTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < nbTilesY; ty++)
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::ComputeFinalTileSize()
{

    uint32_t moduloX = m_ImageWidth  % m_MaxTileSizeX;
    uint32_t moduloY = m_ImageHeight % m_MaxTileSizeY;

    uint32_t maxModuloX = moduloX;
    uint32_t maxModuloY = moduloY;
    uint32_t optMaxTileSizeX = m_MaxTileSizeX;
    uint32_t optMaxTileSizeY = m_MaxTileSizeY;

    for(uint32_t MaxTileSizeX = m_MaxTileSizeX - moduloX; MaxTileSizeX < m_MaxTileSizeX; MaxTileSizeX++)
    {
        moduloX = m_ImageWidth % MaxTileSizeX;
        if(moduloX > maxModuloX)
        {
            maxModuloX = moduloX;
            optMaxTileSizeX = MaxTileSizeX;
        }
    }

    for(uint32_t MaxTileSizeY = m_MaxTileSizeY - moduloY; MaxTileSizeY < m_MaxTileSizeY; MaxTileSizeY++)
    {
        moduloY = m_ImageHeight % m_MaxTileSizeY;
        if(moduloY > maxModuloY)
        {
            maxModuloY = moduloY;
            optMaxTileSizeY = MaxTileSizeY;
        }
    }

    m_MaxTileSizeX = optMaxTileSizeX;
    m_MaxTileSizeY = optMaxTileSizeY;

    m_NumberOfTilesX = m_ImageWidth / m_MaxTileSizeX;
    if(maxModuloX > 0)
    {
        m_NumberOfTilesX++;
    }

    m_NumberOfTilesY = m_ImageHeight / m_MaxTileSizeY;
    if(maxModuloY > 0)
    {
        m_NumberOfTilesY++;
    }
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::WriteGraphIfNecessary(const unsigned int ty,
        const unsigned int tx)
{
    if(m_TileMap.size() > 1)
    {
        // There are other tiles to process, we need to store this one
        // on the local disk of this processor.
        std::stringstream os;
        os << m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
        GraphOperationsType::WriteGraphToDisk(m_Graph, os.str());
        m_Graph->Reset();
    }
}

template< typename TInputImage, typename TOutputGraph>
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::ReadGraphIfNecessary(const unsigned int ty,
        const unsigned int tx)
{
    if(m_TileMap.size() > 1)
    {
        std::stringstream in;
        in << m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
        m_Graph = GraphOperationsType::ReadGraphFromDisk(in.str());
        m_Graph->SetImageWidth(m_ImageWidth);
        m_Graph->SetImageHeight(m_ImageHeight);
        m_Graph->SetNumberOfSpectralBands(m_NumberOfSpectralBands);
        m_Graph->SetOriginX(m_OriginX);
        m_Graph->SetOriginY(m_OriginY);
    }
}

template< typename TInputImage, typename TOutputGraph >
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::CreateOutput()
{
	std::cout << "-------------- CREATE OUTPUT -------------" << std::endl;
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    /************************Write graph on the disk if asked***********************************/
    if(m_SplittedGraph && m_WriteGraph)
    {
    	std::cout << "WRITE SPLITTED GRAPH FOR " << MPIConfig::Instance()->GetMyRank() << std::endl;
        uint32_t tid = 0;
        for(auto& kv : this->m_TileMap)
        {
            uint32_t tx = kv.first % this->m_NumberOfTilesX;
            uint32_t ty = kv.first / this->m_NumberOfTilesX;

            // Load the graph if necessary
            if(this->m_TileMap.size() > 1)
            {
                std::stringstream os;
                os << this->m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
                this->m_Graph = GraphOperationsType::ReadGraphFromDisk(os.str());
                this->m_Graph->SetImageWidth(this->m_ImageWidth);
                this->m_Graph->SetImageHeight(this->m_ImageHeight);
                this->m_Graph->SetNumberOfSpectralBands(this->m_NumberOfSpectralBands);
            }

            std::stringstream os;
            os << this->m_OutputDir << "Graph_" << ty << "_" << tx << ".dat";
            GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());
            tid++;
        }

        mpiConfig->barrier();
    }
    else if(m_WriteGraph && mpiConfig->GetMyRank() == 0)
    {
    	std::cout << "WRITE AGGREGATED GRAPH FOR " << MPIConfig::Instance()->GetMyRank() << std::endl;
        std::stringstream os;
        os << this->m_OutputDir << "Graph_" << 0 << "_" << 0 << ".dat";
        GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());
    }

    /************************Convert graph to image if asked***********************************/
    if(m_WriteLabelImage)
    {
    	//If the graph is splitted, we have to read each temp graph and convert to image
    	if(m_SplittedGraph)
		{
    		for(auto& kv : this->m_TileMap)
			{
				uint32_t tx = kv.first % this->m_NumberOfTilesX;
				uint32_t ty = kv.first / this->m_NumberOfTilesX;
				std::cout << "Convert Graph" << ty << "_" << tx << std::endl;

				std::stringstream os;
				if(this->m_TileMap.size() > 1)
				{
					os << this->m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
					this->m_Graph = GraphOperationsType::ReadGraphFromDisk(os.str());
					this->m_Graph->SetImageWidth(this->m_ImageWidth);
					this->m_Graph->SetImageHeight(this->m_ImageHeight);
					this->m_Graph->SetNumberOfSpectralBands(this->m_NumberOfSpectralBands);
				}

				//Convert to image
				ConvertGraphToImage(ty, tx);
			}
		}
    	else
    	{
    		//Convert to label only if rank == 0 (master proc)
    		if(mpiConfig->GetMyRank() == 0)
    		{
    			//We can directly convert the in-memory graph
    			ConvertGraphToImage(0, 0);
    		}
    	}
    }
}

template< typename TInputImage, typename TOutputGraph >
void
LSImageToGraphScheduler< TInputImage, TOutputGraph>
::FinalGraphAgregation()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    // computing maximum positive value of an int , limitation for sending data with MPI_Send
    int exp = 8*sizeof(int) -1;
    const double dMAX_INT = std::pow(2,exp) - 1;
    const int MAX_INT = static_cast<int>(dMAX_INT);

    std::cout<<"maximum integer supported by system : "<<MAX_INT<<std::endl;

    if(this->m_TileMap.size() > 1)
    {
        // Create the final graph
        this->m_Graph = OutputGraphType::New();
        this->m_Graph->SetImageWidth(this->m_ImageWidth);
        this->m_Graph->SetImageHeight(this->m_ImageHeight);
        this->m_Graph->SetNumberOfSpectralBands(this->m_NumberOfSpectralBands);

        bool hasToProcessDuplicatedNodes = false;

        // Determine the row and col bounds.
        std::unordered_set< uint32_t > rowBounds;
        std::unordered_set< uint32_t > colBounds;

        // Aggregation of its graph
        for(auto& kv : this->m_TileMap)
        {

            uint32_t tx = kv.first % this->m_NumberOfTilesX;
            uint32_t ty = kv.first / this->m_NumberOfTilesX;

            std::cout << " PROCESSOR : " << mpiConfig->GetMyRank() << " , GRAPH AGREGATION FOR " << tx << "/" << ty << std::endl;

            // Check if there are adjacent tiles to this tile which belong to the current processor
            auto neighborTiles = SpatialTools::EightConnectivity(kv.first, this->m_NumberOfTilesX, this->m_NumberOfTilesY);
            for(unsigned short d = 0; d < 8; d++)
            {
                if(neighborTiles[d] > -1 && this->m_TileMap.find(neighborTiles[d]) != this->m_TileMap.end())
                {
                    hasToProcessDuplicatedNodes = true;

                    if(ty > 0)
                    {
                        rowBounds.insert(ty * this->m_MaxTileSizeY);
                        rowBounds.insert(ty * this->m_MaxTileSizeY - 1);
                    }

                    if(ty < this->m_NumberOfTilesY - 1)
                    {
                        rowBounds.insert((ty + 1)* this->m_MaxTileSizeY);
                        rowBounds.insert((ty + 1)* this->m_MaxTileSizeY - 1);
                    }

                    if(tx > 0)
                    {
                        colBounds.insert(tx * this->m_MaxTileSizeX);
                        colBounds.insert(tx * this->m_MaxTileSizeX - 1);
                    }

                    if(tx < this->m_NumberOfTilesX - 1)
                    {
                        colBounds.insert((tx + 1) * this->m_MaxTileSizeX);
                        colBounds.insert((tx + 1) * this->m_MaxTileSizeX - 1);
                    }

                } // end if((neighborTiles[d] > -1 && this->m_TileMap.find(neighborTiles[d]) != this->m_TileMap.end()))

            } // end for(unsigned short d = 0; d < 8; d++)

            std::stringstream os;
            os << this->m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
            auto graph = GraphOperationsType::ReadGraphFromDisk(os.str());
            GraphOperationsType::AggregateGraphs(this->m_Graph, graph);
            graph->Reset();

            if(hasToProcessDuplicatedNodes)
            {
                // Retrieve the nodes on the borders of the adjacent tiles.
                auto borderNodeMap = GraphOperationsType::BuildBorderNodesMapForFinalAggregation(this->m_Graph,
                                                                                                 rowBounds,
                                                                                                 colBounds,
                                                                                                 this->m_ImageWidth);

                // Remove the duplicated nodes
                GraphOperationsType::RemoveDuplicatedNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth);

                // Update the edges
                GraphOperationsType::DetectNewAdjacentNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth, this->m_ImageHeight);
            }

        } // end for(auto& kv : this->m_TileMap)

    } // end if (this->m_TileMap.size() > 1)

    // Synchronisation point
    mpiConfig->barrier();

    std::vector< unsigned long int > sizePerGraph(mpiConfig->GetNbProcs()-1, 0);
    std::vector< char > serializedGraph;

    /* The slave processes have to serialize their graph */
    if(mpiConfig->GetMyRank() == 0)
    {
		//receiving the size of graph for each proc
		for(unsigned int r = 1; r < mpiConfig->GetNbProcs(); r++)
		{
				MPI_Recv(&sizePerGraph[r-1],
							  1,
							  MPI_UNSIGNED_LONG,
							  r,
							  MPI_ANY_TAG,
							  MPI_COMM_WORLD,
							  MPI_STATUS_IGNORE);

				std::cout << " MASTER, size received from proc : " << r << " is " << sizePerGraph[r-1] << std::endl;
		}
    }
    else
    {

        uint32_t tid = 0;
		// Serialize the current graph
		serializedGraph = GraphOperationsType::SerializeGraph(this->m_Graph);

		// sending the size of the graph
		std::cout << " PROCESSOR : " << mpiConfig->GetMyRank() << " , size of the graph : " << serializedGraph.size() << std::endl;
		unsigned long curSize = serializedGraph.size();

		MPI_Send(&curSize,
			 1,
			 MPI_UNSIGNED_LONG,
			 0,
			 0,
			 MPI_COMM_WORLD);

    } // end if (rank !=0 )

    /* Synchronization point */
    mpiConfig->barrier();

    unsigned long int maxNumberOfBytes = serializedGraph.size();

    mpiTools->ComputeMax<unsigned long int>(maxNumberOfBytes, MPI_UNSIGNED_LONG);

    std::cout << "PROC : " << mpiConfig->GetMyRank() << " maxNumberOfBytes " << maxNumberOfBytes <<std::endl;

    if(mpiConfig->GetMyRank() == 0)
    {
        // Simultaneous graph transfer to the master process.
        MPI_Request requests[mpiConfig->GetNbProcs() - 1];
        MPI_Status statuses[mpiConfig->GetNbProcs() - 1];
        //std::vector< std::vector<char> > serializedOtherGraphs(mpiConfig->GetNbProcs() - 1, std::vector<char>(maxNumberOfBytes, char()));
        //std::vector< int > numberOfRecvBytesPerGraph(mpiConfig->GetNbProcs()-1, 0);

	std::vector< std::vector<char> > serializedOtherGraphs;
	serializedOtherGraphs.reserve(mpiConfig->GetNbProcs() - 1);

        for(unsigned int r = 1; r < mpiConfig->GetNbProcs(); r++)
        {

            //serializedOtherGraphs[r-1].assign(maxNumberOfBytes, char());
	    serializedOtherGraphs.push_back(std::vector<char>(sizePerGraph[r-1]));

	    std::cout << "MASTER, receiving data from proc " << r <<std::endl;

	    if(sizePerGraph[r-1] > MAX_INT)
	    {

		// Creating an adapted data_type
		MPI_Datatype severalChar;
		double factor = double(sizePerGraph[r-1]) / dMAX_INT;
		int nbChar = static_cast<int>(ceil(factor));
		MPI_Type_contiguous(nbChar, MPI_CHAR, &severalChar);
		MPI_Type_commit(&severalChar);

		// Computing the count according to the new data_type
		double dCount = double(sizePerGraph[r-1]) / (double)nbChar;
		int recvCount = static_cast<int>(ceil(dCount));

		std::cout << "Receiving " << sizePerGraph[r-1] << " CHAR, factor = " << factor << ", ceil(factor) = "<< ceil(factor)<<", passing by a " << nbChar << " CHAR data_type, counting  "<< recvCount <<" elements. "<<std::endl;

		MPI_Irecv(&(serializedOtherGraphs[r-1][0]),
		              recvCount,
		              severalChar,
		              r,
		              MPI_ANY_TAG,
		              MPI_COMM_WORLD,
		              &(requests[r-1]));
	    }
	    else
	    {
		    std::cout << "Receiving " << sizePerGraph[r-1] << " CHAR, passing by one CHAR data_type."<<std::endl;

		    MPI_Irecv(&(serializedOtherGraphs[r-1][0]),
		              sizePerGraph[r-1],
		              MPI_CHAR,
		              r,
		              MPI_ANY_TAG,
		              MPI_COMM_WORLD,
		              &(requests[r-1]));
	    }

	    // resizing ...
	    std::cout << "MASTER : Resizing for proc " << r << " . Current size = " << serializedOtherGraphs[r-1].size() << ". New size = " << sizePerGraph[r-1] << "...";
	    serializedOtherGraphs[r-1].resize(sizePerGraph[r-1]);
	    std::cout << "OK. " << std::endl;
        }

        if(MPI_Waitall(mpiConfig->GetNbProcs() - 1, requests, statuses) == MPI_ERR_IN_STATUS)
        {
            std::cerr << "Error in processor " << mpiConfig->GetMyRank() << std::endl;
            exit(EXIT_FAILURE);
        }

        // Deserialize and aggregate the graph into mainGraph
        std::cout << "MASTER : Deserialize and aggregate the graph into mainGraph ..."<<std::endl;
        for(unsigned int r = 1; r < mpiConfig->GetNbProcs(); r++)
        {
	    std::cout << "processing proc " << r << " deserialize graph ... ";
            auto otherGraph = GraphOperationsType::DeSerializeGraph(serializedOtherGraphs[r-1]);
	    std::cout << "OK. "<< " clear graph ... ";
            serializedOtherGraphs[r-1].clear();
	    std::cout << "OK. "<< " aggregate graph ... ";
            GraphOperationsType::AggregateGraphs(this->m_Graph, otherGraph);
	    std::cout << "OK. "<< std::endl;

	    std::cout << "test deletion . "<< std::endl;
	    otherGraph->Reset();
        }

        	std::cout << "MASTER : Deserialize and aggregate OK. " << std::endl;

			// Determine the row and col bounds.
        	std::unordered_set< uint32_t > rowBounds;
        	std::unordered_set< uint32_t > colBounds;

			for(uint32_t ty = 0; ty < (unsigned int) this->m_NumberOfTilesY; ty++)
			{
				for(uint32_t tx = 0; tx < (unsigned int) this->m_NumberOfTilesX; tx++)
				{

					if(ty > 0)
					{
						rowBounds.insert(ty * this->m_MaxTileSizeY);
						rowBounds.insert(ty * this->m_MaxTileSizeY - 1);
					}

					if(ty < this->m_NumberOfTilesY - 1)
					{
						rowBounds.insert((ty + 1)* this->m_MaxTileSizeY);
						rowBounds.insert((ty + 1)* this->m_MaxTileSizeY - 1);
					}

					if(tx > 0)
					{
						colBounds.insert(tx * this->m_MaxTileSizeX);
						colBounds.insert(tx * this->m_MaxTileSizeX - 1);
					}

					if(tx < this->m_NumberOfTilesX - 1)
					{
						colBounds.insert((tx + 1) * this->m_MaxTileSizeX);
						colBounds.insert((tx + 1) * this->m_MaxTileSizeX - 1);
					}

				} // end for(uint32_t tx = 0; tx < nbTilesX; tx++)

			} // end for (uint32_t ty = 0; ty < nbTilesY; ty++)


			// Retrieve the nodes on the borders of the adjacent tiles.
			auto borderNodeMap = GraphOperationsType::BuildBorderNodesMapForFinalAggregation(this->m_Graph,
																								 rowBounds,
																								 colBounds,
																								 this->m_ImageWidth);
            // Remove the duplicated nodes
            GraphOperationsType::RemoveDuplicatedNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth);

            // Update the edges
            GraphOperationsType::DetectNewAdjacentNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth, this->m_ImageHeight);

        }
        else
	   {

	    // Passing by MPI_Type_contiguous if the proc exceeds the maximum size of an int
	    if (serializedGraph.size() > MAX_INT)
	    {

			// Creating an adapted data_type
			MPI_Datatype severalChar;
			double factor = double(serializedGraph.size()) / dMAX_INT;
			int nbChar = static_cast<int>(ceil(factor));
			MPI_Type_contiguous(nbChar, MPI_CHAR, &severalChar);
			MPI_Type_commit(&severalChar);

			// Computing the count according to the new data_type
			double dCount = double(serializedGraph.size()) / double(nbChar);
			int currentCount = static_cast<int>(ceil(dCount));

			std::cout << "PROC : " << mpiConfig->GetMyRank() << " , serializedGraph.size() = " << serializedGraph.size() <<" , passing by " << nbChar << "  CHAR, currentCount = "<<currentCount<<std::endl;

			MPI_Send(&serializedGraph[0],
				 currentCount,
				 severalChar,
				 0,
				 0,
				 MPI_COMM_WORLD);
	    }
	    else
	    {

	    	std::cout << "PROC : " << mpiConfig->GetMyRank() << " , serializedGraph.size() = " << serializedGraph.size() <<" , passing by the whole data."<<std::endl;

		    MPI_Send(&serializedGraph[0],
		             serializedGraph.size(),
		             MPI_CHAR,
		             0,
		             0,
		             MPI_COMM_WORLD);
	   }

   }


}

template< typename TInputImage, typename TOutputGraph >
void
LSImageToGraphScheduler<TInputImage, TOutputGraph>
::ConvertGraphToImage(const unsigned int ty, const unsigned int tx)
{
	std::cout << "WRITE LABEL IMAGE FOR " << ty <<"_"<< tx << std::endl;
	if(m_Graph == nullptr)
	{
		std::cout << "Graph " << ty << "_" << tx << " is not set..." << std::endl;
		exit(EXIT_FAILURE);
	}

	using LabelPixelType = unsigned int;
	using LabelImageType = otb::Image< LabelPixelType, 2 >;
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<TOutputGraph, LabelImageType>;
	using WriterType = otb::ImageFileWriter< LabelImageType>;

	/** FOR COLOR IMAGE
	using RGBPixelType = itk::RGBPixel<unsigned char>;
	using RGBImageType = otb::Image<RGBPixelType, 2>;
	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
	*/
	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;

	//Output name
	//std::stringstream os;
	//os << this->m_OutputDir << m_LabelImageName << ty << "_" << tx << ".tif";

	//auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	auto grayWriter = WriterType::New();
	
	/*
	auto fillHoleFilter = FillholeFilterType::New();

	grayWriter->SetFileName(os.str());
	graphToLabelFilter->SetInput(m_Graph);
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	grayWriter->SetInput(fillHoleFilter->GetOutput());
	grayWriter->Update();
	*/

	// Streaming label image writer
  std::stringstream os2;
  os2 << this->m_OutputDir << m_LabelImageName << ty << "_" << tx << "_stream.tif";

  using StreamingGraph2LabelImgFilterType = otb::obia::StreamingGraphToImageFilter<TOutputGraph, LabelImageType>;
  auto filter = StreamingGraph2LabelImgFilterType::New();
  filter->SetInput(m_Graph);

  grayWriter->SetFileName(os2.str());
  grayWriter->SetInput(filter->GetOutput());
  grayWriter->Update();

}

} // end of namespace obia
} // end of namespace otb


#endif
