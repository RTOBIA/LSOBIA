#ifndef otbObiaLSGraphToVectorScheduler_txx
#define otbObiaLSGraphToVectorScheduler_txx

#include "otbObiaLSGraphToVectorScheduler.h"
#include "otbObiaVectorOperations.h"
#include "otbObiaMPITools.h"
namespace otb
{
namespace obia
{

template< typename TInputGraph>
LSGraphToVectorScheduler<TInputGraph>
::LSGraphToVectorScheduler() :
 m_OutputDir(""),
 m_WriteVector(false),
 m_SplittedGraph(false),
 m_WriteGraph(false),
 m_Graph(nullptr),
 m_OutputDS(nullptr),
 m_NumberOfTilesX(0),
 m_NumberOfTilesY(0)
 {
	/**TODO Temporaire, a mettre en paramètre*/
    m_NumberOfTilesX = 2;
    m_NumberOfTilesY = 2;
    m_MaxTileSizeX = 500;
    m_MaxTileSizeY = 500;
    m_ImageHeight = 1000;
    m_ImageWidth = 1000;
//    m_NumberOfTilesX = 1;
//    m_NumberOfTilesY = 1;
//    m_MaxTileSizeX = 1000;
//    m_MaxTileSizeY = 1000;


 }

template< typename TInputGraph>
LSGraphToVectorScheduler<TInputGraph>
::~LSGraphToVectorScheduler()
{

}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::Update()
 {
	//Pre-processing tile frame
    this->PreProcessing();

    //Generate data
    this->GenerateData();

    //Create output
    this->CreateOutput();
 }

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetFileName(const std::string& filename)
{
    m_FileName = filename;
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetOutputDir(const std::string & path)
{
    m_OutputDir = path;
    if( !m_OutputDir.empty() && m_OutputDir[m_OutputDir.size() - 1] != '/')
    {
        m_OutputDir.append("/");
    }
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetMaxTileSizeX(const uint32_t MaxTileSizeX)
{
    m_MaxTileSizeX = MaxTileSizeX;
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetMaxTileSizeY(const uint32_t MaxTileSizeY)
{
    m_MaxTileSizeY = MaxTileSizeY;
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetAvailableMemory(const uint64_t mem)
{
    m_AvailableMemory = mem * 1024 * 1024;
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
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

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::SetInputDirectory(const std::string& inputDir)
{
	m_InputDirectory = inputDir;

	if(!m_InputDirectory.empty())
	{
		if(m_InputDirectory[m_InputDirectory.size()-1] != '/')
		{
			m_InputDirectory.append("/");
		}
	}
}
template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::PreProcessing()
{
	//Padding value is equal to 0
	m_PaddingValue = 0;

	//Compute tile frames
    std::cout << "Number of Tiles Y = " << m_NumberOfTilesY << std::endl;
    std::cout << "Number of Tiles X = " << m_NumberOfTilesX << std::endl;
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
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::ComputeTileFrames()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    uint32_t tid = 0;
    long int rectifiedIndex;
    long int rectifiedPadding;

    //Clear old map
    m_TilesPerProcessor.clear();
    m_TileMap.clear();

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

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::WriteGraphIfNecessary(const unsigned int ty,
					    const unsigned int tx,
						std::string graph_dir,
						std::string prefix_name)
{
    //if(m_TileMap.size() > 1)
    //{
        // There are other tiles to process, we need to store this one
        // on the local disk of this processor.
        std::stringstream os;
        os << graph_dir << prefix_name << "_" << ty << "_" << tx << ".dat";
        GraphOperationsType::WriteGraphToDisk(m_Graph, os.str());
        m_Graph->Reset();
    //}
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::ReadGraphIfNecessary(const unsigned int ty,
					   const unsigned int tx,
					   std::string graph_dir,
					   std::string prefix_name)
{
   // if(m_TileMap.size() > 1)
   // {
        std::stringstream in;
        in << graph_dir << prefix_name << "_" << ty << "_" << tx << ".dat";
        std::cout << "Read graph : "  << in.str() << std::endl;
        m_Graph = GraphOperationsType::ReadGraphFromDisk(in.str());
    //}
}


template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::WriteVectorIfNecessary(const unsigned int ty,
						 const unsigned int tx)
{

	std::cout << "Write If Necessary" << std::endl;

    //if(m_TileMap.size() > 1)
   // {
    	  // There are other tiles to process, we need to store this one
		// on the local disk of this processor.
		std::stringstream os;
		os << m_TemporaryDirectory << "Vector_" << ty << "_" << tx << ".gml";

    	GDALDriver*  poDriverGml = VectorOperations::InitializeGDALDriver("GML");

    	//Create copy
    	GDALDataset* polyGDS = poDriverGml->CreateCopy(os.str().c_str(), &this->m_OutputDS->ogr(), true, nullptr, nullptr, nullptr);

    	//Close polyGDs
    	GDALClose(polyGDS);
   // }
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::ReadVectorIfNecessary(const unsigned int ty,
        const unsigned int tx)
{
    if(m_TileMap.size() > 1)
    {
    }
}

template< typename TInputGraph>
void
LSGraphToVectorScheduler<TInputGraph>
::CreateOutput()
{

}

} // end of namespace obia
} // end of namespace otb


#endif
