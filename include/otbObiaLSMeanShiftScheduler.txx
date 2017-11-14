#ifndef otbObiaLSMeanShiftScheduler_txx
#define otbObiaLSMeanShiftScheduler_txx
#include "otbObiaLSMeanShiftScheduler.h"
#include "otbObiaStreamUtils.h"
#include <unordered_set>

namespace otb
{
namespace obia
{

template< typename TInputImage, typename TLabelPixel>
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::LSMeanShiftScheduler() :
  m_MaxNumberOfIterations(100),
  m_SpatialBandWidth(5),
  m_SpectralRangeBandWidth(15),
  m_Threshold(0.1),
  m_SpectralRangeRamp(0.0),
  m_ModeSearch(true)
{
}

template< typename TInputImage, typename TLabelPixel>
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::~LSMeanShiftScheduler()
{
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::ComputePaddingValue()
{
    /** Stability margin of the MS algorithm + 1 pixel overlap for the CC in order to stitch the tiles. */
    this->m_PaddingValue =  (m_MaxNumberOfIterations * m_SpatialBandWidth) + 1;
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::ComputeExtractionParametersForCC(ProcessingTile& tile,
                                   uint32_t& startX,
                                   uint32_t& startY)
{
    if(tile.m_Ty > 0)
    {
        // A stability margin has been considered at top.
        startY = tile.m_MarginValues[TOP] - 1;

        // Update the y-index of the tile
        tile.m_Frame.SetIndex(1, tile.m_Frame.GetIndex(1) + startY);
      
        // Increment the size of the tile
        tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) + 1);

        // Update the margin value.
        tile.m_MarginValues[TOP] = 1;
    }

    if(tile.m_Tx < this->m_NumberOfTilesX - 1)
    {
        // 1 pixel overlap at right.
        tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) + 1);

        // Update the margin value.
        tile.m_MarginValues[RIGHT] = 1;
    }

    if(tile.m_Ty < this->m_NumberOfTilesY - 1)
    {
        tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) + 1);

        // Update the margin value.
        tile.m_MarginValues[BOTTOM] = 1;
    }

    if(tile.m_Tx > 0)
    {
        // A stability margin has been considered at left
        startX = tile.m_MarginValues[3] - 1;

        // Update the x-index of the tile
        tile.m_Frame.SetIndex(0, tile.m_Frame.GetIndex(0) + startX);

        // 1 pixel overlap at left
        tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) + 1);

        tile.m_MarginValues[LEFT] = 1;
    }    
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::RescaleGraph(ProcessingTile& tile)
{
    auto lambdaOp = [&](NodeType& node)
    {
        node.SetFirstPixelCoords(SpatialTools::TransformPixelCoordsFromTileRefToImgRef(node.GetFirstPixelCoords(),
                                                                                             tile,
                                                                                          this->m_ImageWidth));
        node.m_BoundingBox = SpatialTools::TransformBBoxCoordsFromTileRefToImgRef(node.m_BoundingBox,
                                                                                  tile);

        for(auto& pix : node.m_Attributes.m_ListOfPixels)
        {
            pix = SpatialTools::TransformPixelCoordsFromTileRefToImgRef(pix, tile, this->m_ImageWidth);
        }        
    };

    this->m_Graph->ApplyForEachNode(lambdaOp);

    // Update the tile frame without considering the margins
    tile.m_Frame.SetIndex(0, tile.m_Frame.GetIndex(0) + tile.m_MarginValues[LEFT]);
    tile.m_Frame.SetIndex(1, tile.m_Frame.GetIndex(1) + tile.m_MarginValues[TOP]);
    tile.m_Frame.SetSize(0, tile.m_Frame.GetSize(0) - tile.m_MarginValues[LEFT] - tile.m_MarginValues[RIGHT]);
    tile.m_Frame.SetSize(1, tile.m_Frame.GetSize(1) - tile.m_MarginValues[TOP] - tile.m_MarginValues[BOTTOM]);
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::Segment()
{
  auto mpiConfig = MPIConfig::Instance();
  auto mpiTools = MPITools::Instance();
  
  uint32_t tid = 0;
  
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {

            if(mpiTools->IsMyTurn(tid))
            {

                // Retrieve the current processed tile.
                auto& tile = this->m_TileMap[tid];

                /** Image reader */
                auto imgReader = InputImageFileReaderType::New();
                imgReader->SetFileName(this->m_FileName);

                  /** Extraction of the tile for smoothing step */

                auto tileExtractor = MultiChannelExtractROIFilterType::New();
                tileExtractor->SetStartX(tile.m_Frame.GetIndex(0));
                tileExtractor->SetStartY(tile.m_Frame.GetIndex(1));
                tileExtractor->SetSizeX(tile.m_Frame.GetSize(0));
                tileExtractor->SetSizeY(tile.m_Frame.GetSize(1));

                /** Tile smoothing with the Mean-Shift algorithm */
                auto msFilter = MeanShiftSmoothingImageFilterType::New();
                msFilter->SetSpatialBandwidth(m_SpatialBandWidth);
                msFilter->SetRangeBandwidth(m_SpectralRangeBandWidth);
                msFilter->SetThreshold(m_Threshold);
                msFilter->SetMaxIterationNumber(m_MaxNumberOfIterations);
                msFilter->SetRangeBandwidthRamp(m_SpectralRangeRamp);
                msFilter->SetModeSearch(m_ModeSearch);

                /** Extraction of the stable area within the smooth tile */
                tile.m_Frame.SetSize(0, std::min(this->m_MaxTileSizeX, this->m_ImageWidth - tx * this->m_MaxTileSizeX));
                tile.m_Frame.SetSize(1, std::min(this->m_MaxTileSizeY, this->m_ImageHeight - ty * this->m_MaxTileSizeY));
                uint32_t startX = 0;
                uint32_t startY = 0;
                ComputeExtractionParametersForCC(tile, startX, startY);

                auto stableExtractor = MultiChannelExtractROIFilterType::New();
                tileExtractor->SetStartX(startX);
                tileExtractor->SetStartY(startY);
                tileExtractor->SetSizeX(tile.m_Frame.GetSize(0));
                tileExtractor->SetSizeY(tile.m_Frame.GetSize(1));

                /** Segmentation of the smooth tile with the Connected Component algorithm */
                auto ccFilter = CCFilterType::New();
                std::stringstream expr;
                expr<<"sqrt((p1b1-p2b1)*(p1b1-p2b1)";
                for(unsigned int i=1; i< this->m_NumberOfSpectralBands; i++)
                expr<<"+(p1b"<<i+1<<"-p2b"<<i+1<<")*(p1b"<<i+1<<"-p2b"<<i+1<<")";
                expr<<")"<<"<"<< m_SpectralRangeBandWidth;
                ccFilter->GetFunctor().SetExpression(expr.str());

                /** Extraction of the graph of adjacency from the segmented tile */
                auto graphExtractor = LabelImageToGraphFilterType::New();


                /** Pipeline branching */
                tileExtractor->SetInput(imgReader->GetOutput());
                msFilter->SetInput(tileExtractor->GetOutput());
                stableExtractor->SetInput(msFilter->GetOutput());
                ccFilter->SetInput(stableExtractor->GetOutput());
                ccFilter->Update();
                graphExtractor->SetInput(ccFilter->GetOutput());
                graphExtractor->Update();
                this->m_Graph = graphExtractor->GetOutput();

                /** Tile referential to image referential */
                RescaleGraph(tile);

                // If a processor is in charge of more than one tiles, then this graph has to be stored
                // in the temporary local disk of the execution node.
                this->WriteGraphIfNecessary(ty, tx);

            } //end if(mpiTools->IsMyTurn(tid))

              tid++;

        } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesX; ty++)
    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesX; ty++)

  // Synchronization point
  mpiConfig->barrier();
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::SerializeListOfBorderNodes(const ProcessingTile& tile)
{
    // Get the list of the nodes located at the border of the graph
    auto vecBorderNodes = GraphOperationsType::GetListOfBorderNodes(this->m_Graph, 
                                                                      tile, 
                                                                      this->m_NumberOfTilesX,
                                                                      this->m_NumberOfTilesY,
                                                                      this->m_ImageWidth);

    // Transform this list into a serializable data structure (not so elegant...)
    std::unordered_map< NodeType*, uint32_t > borderNodeMap;
    for(auto& nodePtr : vecBorderNodes)
    {
        borderNodeMap[nodePtr] = 0;
    }
    vecBorderNodes.clear();


    m_SerializedListOfBorderNodes = GraphOperationsType::SerializeStabilityMargin(borderNodeMap, this->m_Graph);

    if(m_MaxNumberOfBytes < m_SerializedListOfBorderNodes.size())
    {
        m_MaxNumberOfBytes = m_SerializedListOfBorderNodes.size();
    }
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::FillSharedBuffer()
{
    auto mpiConfig = MPIConfig::Instance();
      auto mpiTools = MPITools::Instance();

      mpiTools->ComputeMax<unsigned long int>(m_MaxNumberOfBytes, MPI_UNSIGNED_LONG);

      // Allocate the shared buffer
      m_SharedBuffer.clear();
      // Header gives the real number of bytes + the sequence of bytes
      uint64_t maxNumberOfElements = this->m_MaxNumberOfTilesPerProcessor * (IntSize + m_MaxNumberOfBytes);
      m_SharedBuffer.assign(maxNumberOfElements, char());

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
                      m_SerializedListOfBorderNodes = GraphOperationsType::ReadSerializedMarginFromDisk(os.str());
                  }

                  // Move at the right location in the shared buffer.
                  uint64_t offset = ntile * (IntSize + m_MaxNumberOfBytes);

                  // Write the serialized stablity margin in the shared buffer
		  to_stream(m_SharedBuffer,m_SerializedListOfBorderNodes,offset);

                  // We can clean the serialized list of border nodes
                  m_SerializedListOfBorderNodes.clear();

                  ntile++;
            }
            tid++;
        }
    }
    mpiConfig->barrier();
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::AchieveBorderNodeConstruction()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    // Creation of rma window: each processor will have its shared buffer accessible for other processors
    uint64_t maxNumberOfElements = this->m_MaxNumberOfTilesPerProcessor * (IntSize + m_MaxNumberOfBytes);
    MPI_Win win;
    MPI_Win_create(&m_SharedBuffer[0], maxNumberOfElements, CharSize, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    uint32_t tid = 0;
    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {

                // Retrieve the current processed tile.
                auto& tile = this->m_TileMap[tid];

                // Creation of the list of the adjacent serialized list of border nodes. 
                std::vector< std::vector<char> > ajdSerializedListOfBorderNodes;

                // Retrieve the neighbor tiles
                auto neighborTiles = SpatialTools::EightConnectivity(tid, this->m_NumberOfTilesX, this->m_NumberOfTilesY);
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
                        uint64_t offset = pos * (IntSize + m_MaxNumberOfBytes);

                        // Allocate a new serialized stability margin.
                        ajdSerializedListOfBorderNodes.push_back( std::vector<char>(IntSize + m_MaxNumberOfBytes) );

                        // Read rma operation
                        
                        //MPI_Win_fence(0, win);
                        MPI_Win_lock(MPI_LOCK_SHARED, neighRank, 0, win);

                        MPI_Get(&(ajdSerializedListOfBorderNodes[ajdSerializedListOfBorderNodes.size()-1][0]), 
                                IntSize + m_MaxNumberOfBytes,
                                MPI_CHAR,
                                neighRank,
                                offset,
                                IntSize + m_MaxNumberOfBytes,
                                MPI_CHAR,
                                win);

                        MPI_Win_unlock(neighRank, win);
                        //MPI_Win_fence(0, win);

                    } // end if(neighborTiles[n] > -1)

                } // end for(unsigned short n = 0; n < 8; n++)

                // Read the graph if necessary
                this->ReadGraphIfNecessary(ty, tx);

                // Allocate the adjacent subgraph containing the border nodes of the adjacent graph
                std::vector< OutputGraphPointerType > adjSubGraphs(ajdSerializedListOfBorderNodes.size(), nullptr);

                // Deserialize the adjacent subgraphs.
                for(uint32_t i = 0; i < ajdSerializedListOfBorderNodes.size(); i++)
                {
                    // Retrieve the serialized margin
                    std::vector< char > serializedListOfBorderNodes(0);
		    from_stream(ajdSerializedListOfBorderNodes[i],serializedListOfBorderNodes);

                    ajdSerializedListOfBorderNodes[i].clear();
                    adjSubGraphs[i] = GraphOperationsType::DeSerializeGraph(serializedListOfBorderNodes);
                }

                std::unordered_set< uint32_t > rowBounds;
                std::unordered_set< uint32_t > colBounds;

                if(tile.m_Ty > 0)
                {
                    rowBounds.insert(tile.m_Frame.GetIndex(1));
                }

                if(tile.m_Ty < this->m_NumberOfTilesY - 1)
                {
                    rowBounds.insert(tile.m_Frame.GetIndex(1) + tile.m_Frame.GetSize(1) - 1);
                }

                if(tile.m_Tx > 0)
                {
                    colBounds.insert(tile.m_Frame.GetIndex(0));
                }

                if(tile.m_Tx < this->m_NumberOfTilesX - 1)
                {
                    colBounds.insert(tile.m_Frame.GetIndex(0) + tile.m_Frame.GetSize(0) - 1);
                }

                std::unordered_map<CoordValueType, NodeType* > refBorderNodeMap;
                std::vector< std::unordered_map<CoordValueType, NodeType* > > adjBorderNodeMaps(adjSubGraphs.size());
                // Loop over the nodes
                for(auto nodeIt = this->m_Graph->Begin(); nodeIt != this->m_Graph->End(); nodeIt++)
                {
                    // Loop over the border pixels
                    for(auto& pix : nodeIt->m_Attributes.m_ListOfPixels)
                    {
                        if(SpatialTools::IsPixelAtBoundaries(pix,
                                                             rowBounds,
                                                             colBounds,
                                                             this->m_ImageWidth))
                        {
                            refBorderNodeMap[pix] = &(*nodeIt);
                        }

                    } // end for (auto& pix : borderPixels)

                } // end for(auto nodeIt = graph->Begin(); nodeIt != graph->End(); nodeIt++)

                // Deserialize the adjacent subgraphs.
                for(uint32_t i = 0; i < adjSubGraphs.size(); i++)
                {
                    for(auto nodeIt = adjSubGraphs[i]->Begin(); nodeIt != adjSubGraphs[i]->End(); nodeIt++)
                    {
                        // Loop over the border pixels
                        for(auto& pix : nodeIt->m_Attributes.m_ListOfPixels)
                        {
                            if(SpatialTools::IsPixelAtBoundaries(pix,
                                                                 rowBounds,
                                                                 colBounds,
                                                                 this->m_ImageWidth))
                            {
                                adjBorderNodeMaps[i][pix] = &(*nodeIt);
                            }

                        } // end for (auto& pix : borderPixels)

                    }
                }

                // For each pixel in the reference border map, we look for this pixel in all
                // the adjacent border and update the list of internal pixels.
                for(auto& refKV : refBorderNodeMap)
                {
                    // Merge the list of internal pixels.
                    std::unordered_set< CoordValueType > setOfInternalPixels;
                    for(auto& pix : refKV.second->m_Attributes.m_ListOfPixels){ setOfInternalPixels.insert(pix); }

                    for(uint32_t i = 0; i < adjBorderNodeMaps.size(); i++)
                    {

                        auto isHere = adjBorderNodeMaps[i].find(refKV.first);

                        if(isHere != adjBorderNodeMaps[i].end())
                        {
                            auto refNode = refKV.second;
                            auto complNode = isHere->second;

                            for(auto& pix : complNode->m_Attributes.m_ListOfPixels)
                            {
                                if(setOfInternalPixels.find(pix) == setOfInternalPixels.end())
                                {
                                    refNode->m_Attributes.m_ListOfPixels.push_back(pix);
                                    setOfInternalPixels.insert(pix);
                                }
                            }

                            // Fusion of the contour
                            refNode->m_Contour.MergeWith(complNode->m_Contour, this->m_ImageWidth, this->m_ImageHeight);

                        }

                    }

                }

                // If a processor is in charge of more than one tiles, then this graph has to be stored
                // in the temporary local disk of the execution node.
                this->WriteGraphIfNecessary(ty, tx);

            }

            tid++;
        }
    }


    mpiConfig->barrier();

    // Can release the rma window
    MPI_Win_free(&win);

    m_SharedBuffer.clear();
}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::PostProcessing()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();
  
    uint32_t tid = 0;
    m_MaxNumberOfBytes = 0;

    for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)
    {
        for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)
        {
            if(mpiTools->IsMyTurn(tid))
            {

                // Retrieve the current processed tile.
                auto& tile = this->m_TileMap[tid];

                // Read the graph if necessary
                this->ReadGraphIfNecessary(ty, tx);

                // Extract the border nodes and serialize them
                SerializeListOfBorderNodes(tile);

                if(this->m_TileMap.size() > 1)
                {
                    std::stringstream os;
                    os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
                    GraphOperationsType::WriteSerializedMarginToDisk(m_SerializedListOfBorderNodes, os.str());
                }
            } // end if(mpiTools->IsMyTurn(tid))

            tid++;

        } // end for(uint32_t tx = 0; tx < this->m_NumberOfTilesX; tx++)

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)

    // Synchronization point
    mpiConfig->barrier();

    // Build the remote window.
    FillSharedBuffer();

    // Build a border pixel map: Map: pixel -> list of nodes
    AchieveBorderNodeConstruction();

}

template< typename TInputImage, typename TLabelPixel>
void
LSMeanShiftScheduler<TInputImage, TLabelPixel>
::GenerateData()
{
    /** Mean-Shift -> Connected Component -> Graph Extraction */
    Segment();

    /** Achieve the construction of the nodes located at the borders of the graphs */
    PostProcessing();
}


} // end of namespace obia
} // end of namespace otb

#endif
