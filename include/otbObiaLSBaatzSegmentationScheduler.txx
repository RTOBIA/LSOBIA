#ifndef otbObiaLSBaatzSegmentationScheduler_txx
#define otbObiaLSBaatzSegmentationScheduler_txx
#include "otbObiaLSBaatzSegmentationScheduler.h"
#include "otbObiaGraphToLabelImageFilter.h"
#include "otbObiaStreamUtils.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImage.h"
#include "otbImageFileWriter.h"
#include <unordered_set>


namespace otb
{
namespace obia
{

template< class TInputImage >
LSBaatzSegmentationScheduler<TInputImage>
::LSBaatzSegmentationScheduler() :
m_MaxNumberOfIterations(75),
m_CurrentNumberOfIterations(0),
m_Threshold(1600),
m_Decreasing(0.0),
m_SpectralWeight(0.5),
m_ShapeWeight(0.5),
m_StartingNumberOfIterations(1),
m_AggregateGraphs(false),
m_MaxNumberOfBytes(0),
m_PartialNumberOfIterations(1)
{
    m_BandWeights.clear();
}

template< class TInputImage >
LSBaatzSegmentationScheduler<TInputImage>
::~LSBaatzSegmentationScheduler()
{
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::ComputePaddingValue()
{
    /** Stability margin of the large scale generic region merging algorithm */
    this->m_PaddingValue =  pow(2, m_StartingNumberOfIterations + 1) - 2;
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::NoTilingExecution()
{
	std::cout << "NO TILING EXECUTION WITH " << m_MaxNumberOfIterations << std::endl;
    auto mpiConfig = MPIConfig::Instance();

    // Only the master node segments the whole input image.
    if(mpiConfig->GetMyRank() == 0)
    {
        // Read the input image
        auto imgReader = InputImageReaderType::New();
        imgReader->SetFileName(this->m_FileName);

        // Creation of the initial baatz graph
        auto imgToBaatzFilter = ImageToBaatzGraphFilterType::New();

        // Segmentation filter
        auto baatzFilter = CreateFilter();
        baatzFilter->SetMaxNumberOfIterations(m_MaxNumberOfIterations);

        // Pipeline branching
        imgToBaatzFilter->SetInput(imgReader->GetOutput());

	imgToBaatzFilter->SetProcessNoData(this->m_ProcessNoData);
	imgToBaatzFilter->SetNoDataValue(this->m_NoDataValue);
        baatzFilter->SetInput(imgToBaatzFilter->GetOutput());
        baatzFilter->Update();

        this->m_Graph = baatzFilter->GetOutput();

    }

    this->m_SplittedGraph = false;
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::TilingExecution()
{
	std::cout << "TILING EXECUTION WITH " << m_StartingNumberOfIterations << std::endl;

    // First partial segmentation
    auto segState = FirstPartialSegmentation();

    m_CurrentNumberOfIterations = m_StartingNumberOfIterations;

    if(m_CurrentNumberOfIterations < m_MaxNumberOfIterations)
    {
        switch(segState)
        {
            case PARTIAL_SEG:
              std::cout << "PARTIAL_SEG" << std::endl;
              PartialSegmentation();
            break;
            case SEG_OVER_NO_AGGR:
              std::cout << "SEG_OVER_NO_AGGR" << std::endl;
              this->m_SplittedGraph = true;
            break;
            case AGGR_AND_SEG:
              std::cout << "AGGR_AND_SEG" << std::endl;
              FinalGraphAgregation();
              AchieveSegmentation();
              this->m_SplittedGraph = false;
            break;
            case AGGR_NO_SEG:
              std::cout << "AGGR_NO_SEG" << std::endl;
              FinalGraphAgregation();
              this->m_SplittedGraph = false;
            break;
            case NO_AGGR_AND_SEG:
            std::cout << "NO_AGGR_AND_SEG" << std::endl;
            NoAggregationAndPartialSegmentation();
            this->m_SplittedGraph = true;
            break;
            case SEG_STATE_UNDEF:
            std::cout << "SEG_STATE_UNDEF, should stop" << std::endl;
            break;
            default:
                std::cout << "DEFAULT state, should stop" << std::endl;

        }
    }
}

template< class TInputImage >
itk::SmartPointer<typename LSBaatzSegmentationScheduler<TInputImage>::BaatzSegmentationFilterType> LSBaatzSegmentationScheduler<TInputImage>
::CreateFilter()
{
    //Create the 3 classes
    auto mergingFunc = new BaatzMergingCost<float, GraphType>(); //::New();
    mergingFunc->SetSpectralWeight(this->m_SpectralWeight);
    mergingFunc->SetShapeWeight(this->m_ShapeWeight);
    mergingFunc->SetBandWeights(this->m_BandWeights);
    mergingFunc->SetThreshold(this->m_Threshold);

    auto heuristicFunc = new BaatzHeuristic<GraphType>();//::New();
    heuristicFunc->SetGraph(this->GetGraph());
    heuristicFunc->SetThreshold(this->m_Threshold);

    auto updateFunc = new BaatzUpdateAttribute<GraphType>();//::New();

    // Segmentation filter
    auto baatzFilter = BaatzSegmentationFilterType::New();
    baatzFilter->SetMaxNumberOfIterations(this->m_MaxNumberOfIterations);
    baatzFilter->SetMergingCostFunc(mergingFunc);
    baatzFilter->SetHeuristicFunc(heuristicFunc);
    baatzFilter->SetUpdateAttributeFunc(updateFunc);

    return baatzFilter;
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::GenerateData()
{
    /** Reset the current number of iterations applied */
    m_CurrentNumberOfIterations = 0;


    // If the number of processors is greater than the number of processors
    // which have actually a charge than returns an error indicating the real
    // number of processes to launch with the MPI command.
    // If there is just one processor with one tile then the master node segments
    // the image entirely.
    auto mpiConfig = otb::MPIConfig::Instance();

    if(this->m_TilesPerProcessor.size() == 1 && this->m_TilesPerProcessor.begin()->second.size() == 1)
    {
        // Only the master node segments the entire input image.
        NoTilingExecution();
    }
    else if(this->m_TilesPerProcessor.size() < mpiConfig->GetNbProcs())
    {
        std::cerr << "Error: You have to launch your MPI program with " << this->m_TilesPerProcessor.size() << " processor(s)." << std::endl;
        mpiConfig->abort(1);
    }
    else
    {
        TilingExecution();
    }


}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::RescaleGraph(ProcessingTile& tile)
{
    auto lambdaOp = [&](NodeType& node)
    {
        node.SetFirstPixelCoords(SpatialTools::TransformPixelCoordsFromTileRefToImgRef(node.GetFirstPixelCoords(),
                                                                                             tile,
                                                                                          this->m_ImageWidth));
        node.m_BoundingBox = SpatialTools::TransformBBoxCoordsFromTileRefToImgRef(node.m_BoundingBox,
                                                                                  tile);
    };

    this->m_Graph->ApplyForEachNode(lambdaOp);

    // Update the tile frame without considering the margins
    tile.m_Frame.SetIndex(0, tile.m_Frame.GetIndex(0) + tile.m_MarginValues[LEFT]);
    tile.m_Frame.SetIndex(1, tile.m_Frame.GetIndex(1) + tile.m_MarginValues[TOP]);
    tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) - tile.m_MarginValues[LEFT] - tile.m_MarginValues[RIGHT]);
    tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) - tile.m_MarginValues[TOP] - tile.m_MarginValues[BOTTOM]);
}

template< class TInputImage >
typename LSBaatzSegmentationScheduler<TInputImage>::SegState
LSBaatzSegmentationScheduler<TInputImage>
::FirstPartialSegmentation()
{
	std::cout <<"----------- First Partial Segmentation with " << m_StartingNumberOfIterations << std::endl;
	auto mpiConfig = MPIConfig::Instance();
	auto mpiTools = MPITools::Instance();

	uint32_t tid = 0;
	int localFusion = 0;
	unsigned long int accumulatedMemory = 0;

	for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
	{
		for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
		{
			if(mpiTools->IsMyTurn(tid))
			{
				std::cout << "Ty = " << ty  << " Tx = " << tx << std::endl;
				// Retrieve the tile by reference since it will be modified.
				auto& tile = this->m_TileMap[tid];

        // Read the input image
        auto imgReader = InputImageReaderType::New();
        imgReader->SetFileName(this->m_FileName);

        typename InputImageType::PointType origin = imgReader->GetOutput()->GetOrigin();

				// Extraction of the tile
				auto tileExtractor = MultiChannelExtractROIFilterType::New();
				tileExtractor->SetStartX(tile.m_Frame.GetIndex(0));
				tileExtractor->SetStartY(tile.m_Frame.GetIndex(1));
				tileExtractor->SetSizeX(tile.m_Frame.GetSize(0));
				tileExtractor->SetSizeY(tile.m_Frame.GetSize(1));
        tileExtractor->SetInput(imgReader->GetOutput());
        tileExtractor->Update();
        typename InputImageType::Pointer img = tileExtractor->GetOutput();

        // Clear the reader
        imgReader=ITK_NULLPTR;
        tileExtractor=ITK_NULLPTR;

				// Creation of the initial baatz graph
				auto imgToBaatzFilter = ImageToBaatzGraphFilterType::New();

        // Pipeline branching
				imgToBaatzFilter->SetInput(img);
				imgToBaatzFilter->SetProcessNoData(this->m_ProcessNoData);
				imgToBaatzFilter->SetNoDataValue(this->m_NoDataValue);
        imgToBaatzFilter->Update();

				//Update image origin in graph
				this->m_Graph = imgToBaatzFilter->GetOutput();

        // Now clear the tileExtractor and baatzFilter
        tileExtractor=ITK_NULLPTR;
        imgToBaatzFilter=ITK_NULLPTR;

				// Segmentation filter
				// Baatz & Shäpe segmentation
				auto baatzFilter = CreateFilter();
				baatzFilter->SetMaxNumberOfIterations(this->m_StartingNumberOfIterations);
				baatzFilter->SetInput(this->m_Graph);
				baatzFilter->Update();

				// Determine if the segmentation is over
				localFusion += (baatzFilter->GetMergingOver()) ? 0 : 1;

				//Update image origin in graph
				this->m_Graph = baatzFilter->GetOutput();

        // Now clear the baatzFilter
        baatzFilter=ITK_NULLPTR;

				// tile referential -> image referential
				RescaleGraph(tile);

				// Now we are in the image referential
				this->m_Graph->SetImageWidth(this->m_ImageWidth);
				this->m_Graph->SetImageHeight(this->m_ImageHeight);
				this->m_Graph->SetNumberOfSpectralBands(this->m_NumberOfSpectralBands);
				this->m_Graph->SetOriginX(origin[0]);
				this->m_Graph->SetOriginY(origin[1]);

				// Remove the unstable segments
				GraphOperationsType::RemoveUnstableNodes(this->m_Graph,
												         tile,
												         this->m_ImageWidth);

				// Get the memory size of the graph.
				accumulatedMemory += this->m_Graph->GetMemorySize();
				
				// Write the graph if necessary
				this->WriteGraphIfNecessary(ty, tx);

			} // end if(mpiTools->IsMyTurn(tid))

			tid++;

		} // end for(uint32_t tx = 0; tx < nbTilesX; tx++)

	} // end for(uint32_t ty = 0; ty < nbTilesY; ty++)

	/* Synchronization point */
	mpiConfig->barrier();

	/* Compute the accumulated memory */
	mpiTools->Accumulate(accumulatedMemory,  MPI_UNSIGNED_LONG);

	/* Determine if segmentation is globally over */
	mpiTools->Accumulate(localFusion,  MPI_INT);

	if(accumulatedMemory > this->m_AvailableMemory && localFusion > 0)
	{
		return PARTIAL_SEG;
	}
	else if(accumulatedMemory > this->m_AvailableMemory && localFusion < 1)
	{
		return SEG_OVER_NO_AGGR;
	}
	else if(accumulatedMemory < this->m_AvailableMemory && localFusion > 0)
	{
		if(m_AggregateGraphs)
			return AGGR_AND_SEG;
		else
			return NO_AGGR_AND_SEG;
	}
	else if(accumulatedMemory < this->m_AvailableMemory && localFusion < 1)
	{
		if(m_AggregateGraphs)
			return AGGR_NO_SEG;
		else
			return SEG_OVER_NO_AGGR;
	}else
	{
		return SEG_STATE_UNDEF;
	}

}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::PartialSegmentation()
{
    int fusionSum = 1;
    unsigned long int accumulatedMemory = this->m_AvailableMemory + 1;
    unsigned long int prevAccumulatedMemory;
    bool keepOn = true;
    bool firstIt = true;

    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    //initialize margins
    ExtractStabilityMargins();

    while( m_CurrentNumberOfIterations <= m_MaxNumberOfIterations &&
           accumulatedMemory > this->m_AvailableMemory &&
           fusionSum > 0 && keepOn)
    {
        std::cout << "Current iteration " << m_CurrentNumberOfIterations << "/" << m_MaxNumberOfIterations << std::endl;
        std::cout << "Memory : " << accumulatedMemory << "/" << this->m_AvailableMemory << std::endl;
        std::cout << "Local sum : " << fusionSum << std::endl;
        //TODO : ExtractStabilityMargins();

	prevAccumulatedMemory = accumulatedMemory;

        AggregateStabilityMargins();
        RunPartialSegmentation(accumulatedMemory, fusionSum);
        m_CurrentNumberOfIterations += m_PartialNumberOfIterations;

        //Extract for next iterations
        ExtractStabilityMargins();

	//Test decreasing of accumulated Memory : not for first iteration
	if (!firstIt && m_Decreasing)
	{
		float decreasing = (float)(prevAccumulatedMemory - accumulatedMemory)*100.0 / (float)prevAccumulatedMemory;
		std::cout << "prevAccumulatedMemory : " << prevAccumulatedMemory << std::endl;
		std::cout << "accumulatedMemory : " << accumulatedMemory << std::endl;
		std::cout << "decreasing since last iteration : " << decreasing << std::endl;
		if (decreasing < m_Decreasing)
		{
			keepOn = false;
			std::cout << "stopping partial segmentation"<<std::endl;
		}
	}

	firstIt = false;	
	
    } // end while( accumulatedMemory > this->m_AvailableMemory && fusionSum > 0)

    if(accumulatedMemory < this->m_AvailableMemory)
    {
        std::cout << "Aggregation (accumulated memory = " << accumulatedMemory << " and available = "
        		  << this->m_AvailableMemory << ")" << std::endl;
        FinalGraphAgregation();
        this->m_SplittedGraph = false;

        if(m_CurrentNumberOfIterations < m_MaxNumberOfIterations && fusionSum > 0)
        {
            std::cout << "Ending segmentation " << std::endl;
            AchieveSegmentation();
        }
    }
    else
    {
        // For now, we do not write the label image in this situation.
        //this->m_WriteLabelImage = false;

        // Since we cannot aggregate the graph is splitted
        this->m_SplittedGraph = true;
    }

}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::ExtractStabilityMargins()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    // Compute the number of adjacency layers to extract
    const uint32_t nbAdjacencyLayers = pow(2, m_PartialNumberOfIterations + 1) - 2;

    // the local serialized margin
    m_MaxNumberOfBytes = 0;

    uint32_t tid = 0;
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {

                // Load the graph if necessary
                this->ReadGraphIfNecessary(ty, tx);

                // Retrieve the tile
                auto tile = this->m_TileMap[tid];

                auto subGraphMap = GraphOperationsType::ExtractStabilityMargin(this->m_Graph,
                                                                                  nbAdjacencyLayers,
                                                                                  tile,
                                                                                  this->m_NumberOfTilesX,
                                                                                  this->m_NumberOfTilesY,
                                                                                  this->m_ImageWidth
                                                                                  /*this->m_ImageHeight*/);

                m_SerializedStabilityMargin = GraphOperationsType::SerializeStabilityMargin(subGraphMap,
                                                                                              this->m_Graph);

                if(m_MaxNumberOfBytes < m_SerializedStabilityMargin.size())
                {
                    m_MaxNumberOfBytes = m_SerializedStabilityMargin.size();
                }


                if(this->m_TileMap.size() > 1)
                {
                    std::stringstream os;
                    os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
                    GraphOperationsType::WriteSerializedObjectToDisk(m_SerializedStabilityMargin, os.str());
                    m_SerializedStabilityMargin.clear();
                }

            } // end if(mpiTools->IsMyTurn(tid))

            tid++;

        } // end for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)

    mpiConfig->barrier();

    mpiTools->ComputeMax<unsigned long int>(m_MaxNumberOfBytes, MPI_UNSIGNED_LONG);
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::AggregateStabilityMargins()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();


    // Create the shared buffer which will be accessible by other processes.
    uint64_t blockSize = sizeof(size_t) + m_MaxNumberOfBytes;
    uint64_t maxNumberOfElements = this->m_MaxNumberOfTilesPerProcessor * blockSize;
    std::vector< char > sharedBuffer(maxNumberOfElements);

    uint32_t tid = 0;
    uint64_t ntile = 0;
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {

                if(this->m_TileMap.size() > 1)
                {

                    std::stringstream os;
                    os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
                    m_SerializedStabilityMargin = GraphOperationsType::ReadSerializedObjectFromDisk(os.str());
                }

                // Move at the right location in the shared buffer.
                uint64_t offset = ntile * blockSize;

                // Write the serialized stablity margin in the shared buffer
                to_stream(sharedBuffer,m_SerializedStabilityMargin,offset);
                uint64_t numNodes;
                // Write the number of bytes in the serialized margin.
                from_stream(m_SerializedStabilityMargin,numNodes);

                // Can release this serialized stability margin

                std::vector<char>().swap(m_SerializedStabilityMargin);
                    //m_SerializedStabilityMargin.clear();
                    // m_SerializedStabilityMargin.shrink_to_fit();


                // Increment the number of tiles processed
                ntile++;

            } // end if(mpiTools->IsMyTurn(tid))

            tid++;

        } // end for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)


    // Creation of rma window: each processor will have its shared buffer accessible for other processors
    MPI_Win win;
    int success = MPI_Win_create(&sharedBuffer[0], maxNumberOfElements, CharSize, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    assert(success==0);


    tid = 0;
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {
                // Creation of the list of serialized stability margins per processor
                std::vector< std::vector<char> > otherSerializedMargins;

                // Retrieve the neighbor tiles
                auto neighborTiles = SpatialTools::EightConnectivity(tid, this->m_NumberOfTilesX, this->m_NumberOfTilesY);

                MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), win);
                for(unsigned short n = 0; n < 8; n++)
                {
                    if(neighborTiles[n] > -1)
                    {
                        int neighRank = mpiTools->GetProcessorRankFromTileId(neighborTiles[n]);

                        // Retrieve the position of the neigh tile in the shared buffer of the
                        // processor neighRank in order to compute the offset of displacement.
                        uint32_t pos = 0;
                        for(auto& t : this->m_TilesPerProcessor[neighRank])
                        {
                            if(t == neighborTiles[n])
                            {
                                break;
                            }
                            else
                            {
                                pos++;
                            }
                        }

                        // Compute the offset of displacement.
                        uint64_t offset = pos * blockSize;

                        // Allocate a new serialized stability margin.
                        otherSerializedMargins.push_back( std::vector<char>(blockSize) );


                        success = MPI_Get(&(otherSerializedMargins.back()[0]),
                        		blockSize,
                                MPI_CHAR,
                                neighRank,
                                offset,
                                blockSize,
                                MPI_CHAR,
                                win);
                        assert(success==0);





                    } // end if(neighborTiles[n] > -1)

                } // end for(unsigned short n = 0; n < 8; n++)

                 MPI_Win_fence(MPI_MODE_NOSUCCEED,win);

                this->ReadGraphIfNecessary(ty, tx);

                // Agregate the stability margins to the graph
                for(uint32_t i = 0; i < otherSerializedMargins.size(); i++)
                {
   		    std::vector<char> otherSerializedMargin;
                    from_stream(otherSerializedMargins[i],otherSerializedMargin);


                    std::vector<char>().swap(otherSerializedMargins[i]);
                    // otherSerializedMargins[i].clear();
                    // otherSerializedMargins[i].shrink_to_fit();

		    uint64_t numNodes;
		    from_stream(otherSerializedMargin,numNodes);

                    // Deserialize the graph

		    // TODO: This is the call that fails with the duplicated node error
                    auto subGraph = GraphOperationsType::DeSerializeGraph(otherSerializedMargin);

                    GraphOperationsType::AggregateGraphs(this->m_Graph, subGraph);

                    //Reset subgraph
            		subGraph->Reset();
                }

                // Retrieve the tile
                auto tile = this->m_TileMap[tid];

                // Remove duplicated nodes
                auto borderNodeMap = GraphOperationsType::BuildBorderNodesMap(this->m_Graph,
                                                                              tile,
                                                                              this->m_NumberOfTilesX,
                                                                              this->m_NumberOfTilesY,
                                                                              this->m_ImageWidth);

                GraphOperationsType::RemoveDuplicatedNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth);


                // Update edges
                GraphOperationsType::DetectNewAdjacentNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth, this->m_ImageHeight);


                this->WriteGraphIfNecessary(ty, tx);


            } // end if(mpiTools->IsMyTurn(tid))

            tid++;

        } // end for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)


    // Can release the rma window
    MPI_Win_free(&win);
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::RunPartialSegmentation(unsigned long int& accumulatedMemory, int& fusionSum)
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    accumulatedMemory = 0;
    fusionSum = 0;
    uint32_t tid = 0;
    std::cout << "Run Partial Segmentation" << std::endl;
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {

                // Load the graph if necessary
                this->ReadGraphIfNecessary(ty, tx);

                // Retrieve the tile
                auto tile = this->m_TileMap[tid];


                // Segmentation filter
                auto baatzFilter = CreateFilter();
                baatzFilter->SetMaxNumberOfIterations(std::min(m_PartialNumberOfIterations, m_MaxNumberOfIterations - m_CurrentNumberOfIterations + 1));

                // Pipeline branching
                baatzFilter->SetInput(this->m_Graph);
                baatzFilter->Update();
                this->m_Graph = baatzFilter->GetOutput();

                // Determine if the segmentation is over
                fusionSum += (baatzFilter->GetMergingOver()) ? 0 : 1;

                // Remove the unstable segments
                GraphOperationsType::RemoveUnstableNodes(this->m_Graph,
                                                         tile,
                                                         this->m_ImageWidth);

                accumulatedMemory += this->m_Graph->GetMemorySize();

                this->WriteGraphIfNecessary(ty, tx);

            } // end if(mpiTools->IsMyTurn(tid))

            tid++;

        } // end for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)

    mpiConfig->barrier();

    /* Compute the accumulated memory */
    mpiTools->Accumulate(accumulatedMemory,  MPI_UNSIGNED_LONG);

    /* Determine if segmentation is globally over */
    mpiTools->Accumulate(fusionSum,  MPI_INT);
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::FinalGraphAgregation()
{
	SuperClass::FinalGraphAgregation();
}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::NoAggregationAndPartialSegmentation()
{
	std::cout << "-------------NO AGGR AND PARTIAL SEG--------------" << std::endl;
    int fusionSum = 1;
    unsigned long int accumulatedMemory = this->m_AvailableMemory + 1;

    //Initialize stability margins
    ExtractStabilityMargins();

    /**TODO : Add condition when localSum < 1*/
    while(m_CurrentNumberOfIterations <= m_MaxNumberOfIterations &&
          fusionSum > 0)
    {
    	std::cout << "Current iteration = " << m_CurrentNumberOfIterations << std::endl;
       // ExtractStabilityMargins();

        AggregateStabilityMargins();

        RunPartialSegmentation(accumulatedMemory, fusionSum);

        m_CurrentNumberOfIterations += m_PartialNumberOfIterations;

        //Extract stability margins for next iteration
        ExtractStabilityMargins();
    }

}

template< class TInputImage >
void
LSBaatzSegmentationScheduler<TInputImage>
::AchieveSegmentation()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();
    std::cout<<"----------------------"<<std::endl;
    std::cout<<"---ACHIEVING SEGMENTATION -"<<std::endl;
    std::cout<<"----------------------"<<std::endl;

    if(mpiConfig->GetMyRank() == 0 && m_AggregateGraphs)
    {
        std::cout << "Baatz with " << m_MaxNumberOfIterations + 1 - m_CurrentNumberOfIterations << std::endl;
        auto baatzFilter = CreateFilter();
        //auto baatzFilter = BaatzSegmentationFilterType::New();
        baatzFilter->SetInput(this->m_Graph);
        baatzFilter->SetMaxNumberOfIterations(m_MaxNumberOfIterations + 1 - m_CurrentNumberOfIterations);
        /*baatzFilter->SetThreshold(m_Threshold);
        baatzFilter->SetSpectralWeight(m_SpectralWeight);
        baatzFilter->SetShapeWeight(m_ShapeWeight);
        baatzFilter->SetBandWeights(m_BandWeights);*/

        baatzFilter->Update();
        this->m_Graph = baatzFilter->GetOutput();
        this->m_Graph->SetProjectionRef(this->m_ProjectionRef);
        this->m_Graph->SetOriginX(this->m_OriginX);
        this->m_Graph->SetOriginY(this->m_OriginY);
    }
}



} // end of namespace obia
} // end of namespace otb

#endif
