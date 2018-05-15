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
