#ifndef otbObiaLSPolygonizeScheduler_txx
#define otbObiaLSPolygonizeScheduler_txx

#include "otbObiaLSPolygonizeScheduler.h"
#include "otbObiaMPITools.h"
#include <unordered_set>

#include "itkRGBPixel.h"
#include "otbImageFileWriter.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbObiaGraphToLabelImageFilter.h"

namespace otb
{
namespace obia
{

template< class TInputGraph, class TSimplifyFunc >
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::LSPolygonizeScheduler() :
m_AggregateGraphs(false)
{
}

template< class TInputGraph, class TSimplifyFunc >
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::~LSPolygonizeScheduler()
{
	std::cout << "Destruction" << std::endl;
	if(m_SimplifyFunc != nullptr)
	{
		delete m_SimplifyFunc;
	}

	if(this->m_Graph != nullptr)
	{
		this->m_Graph->Reset();
	}

	if(this->m_OutputDS != nullptr)
	{
		this->m_OutputDS->Clear();
	}
}


template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::GenerateData()
{


    // If the number of processors is greater than the number of processors
    // which have actually a charge than returns an error indicating the real
    // number of processes to launch with the MPI command.
    // If there is just one processor with one tile then the master node segments
    // the image entirely.
    auto mpiConfig = otb::MPIConfig::Instance();

    if(this->m_TilesPerProcessor.size() == 1 && this->m_TilesPerProcessor.begin()->second.size() == 1)
    {
    	//Set false for tile processing
    	m_IsTileProcessing = false;

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
    	//Set true for tile processing
    	m_IsTileProcessing = true;

        TilingExecution();
    }
}


template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::NoTilingExecution()
{
    auto mpiConfig = MPIConfig::Instance();
	std::cout << "NO TILING EXECUTION WITH " << mpiConfig->GetMyRank() << std::endl;

    // Only the master node segments the whole input image.
    if(mpiConfig->GetMyRank() == 0)
    {
    	//Temp
    	this->ReadGraphIfNecessary(0, 0, this->m_InputDirectory, "Graph");

    	//Create and update filter
    	RunFilters();

    	//Write vector
    	this->WriteVectorIfNecessary(0, 0);
    }
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::TilingExecution()
{
	std::cout << "TILING EXECUTION " << std::endl;
    auto mpiTools = MPITools::Instance();

    /**************************************************************************/
	/***********************DEBUG : LABEL IMAGE IN RGB ************************/
    /**************************************************************************/
    for(auto& kv : this->m_TileMap)
	{
    	uint32_t tx = kv.first % this->m_NumberOfTilesX;
		uint32_t ty = kv.first / this->m_NumberOfTilesX;

		//Current tile
		m_CurrentTile = kv.second;

    	// Load the graph if necessary
		this->ReadGraphIfNecessary(ty, tx, this->m_InputDirectory, this->m_GraphPrefixName);

		std::stringstream os;
		os.clear();
		os << this->m_OutputDir << "BEFORE_MARGIN_PROC_" << MPIConfig::Instance()->GetMyRank()
		   << "_Tile_" << ty << "_" << tx << ".tif";

		ConvertGraphToImage(this->m_Graph, os.str());
	}

    /**************************************************************************/
	/**************************END DEBUG***********************************************/
    /**************************************************************************/

	//Extract stability margin
	ExtractStabilityMargins();

	//Aggregate graph
	AggregateStabilityMargins();


	/**************************************************************************/
	/***********************DEBUG : LABEL IMAGE IN RGB ************************/
	/**************************************************************************/
	for(auto& kv : this->m_TileMap)
	{
		uint32_t tx = kv.first % this->m_NumberOfTilesX;
		uint32_t ty = kv.first / this->m_NumberOfTilesX;

		// Load the graph if necessary
		this->ReadGraphIfNecessary(ty, tx, this->m_TemporaryDirectory, "Graph_With_Marge");

		std::stringstream os;
		os.clear();
		os << this->m_OutputDir << "AFTER_MARGIN_PROC_" << MPIConfig::Instance()->GetMyRank()
		   << "_Tile_" << ty << "_" << tx << ".tif";

		ConvertGraphToImage(this->m_Graph, os.str());

	}

	/**************************************************************************/
	/**************************END DEBUG***********************************************/
	/**************************************************************************/
	//Run filter on each tile
    for(auto& kv : this->m_TileMap)
    {
        uint32_t tx = kv.first % this->m_NumberOfTilesX;
        uint32_t ty = kv.first / this->m_NumberOfTilesX;

    	std::cout << "Tile " << tx << "/" << ty << std::endl;
		// Load the graph if necessary
		this->ReadGraphIfNecessary(ty, tx, this->m_TemporaryDirectory, "Graph_With_Marge");
		//this->ReadGraphIfNecessary(ty, tx, this->m_InputDirectory, this->m_GraphPrefixName);
		//Run filters
		RunFilters();

		//Remove outside polygons
		//RemovePolygonsOutsideTile(kv.second);

		if(this->m_WriteVector)
		{
			WriteFeatures(kv.second);
		}
    }

    //Wait all
    MPIConfig::Instance()->barrier();

    std::cout << "All filters over" << std::endl;
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::RunFilters()
{
	std::cout << "------------------ RUN FILTERS -----------------" << std::endl;
	//initialize filter (graphToVector)
	//Create filter to vectorize
	auto graphToVectorFilter = GraphToVectorFilterType::New();

	std::cout << "Number of nodes in the graph : " << this->m_Graph->GetNumberOfNodes() << std::endl;
	graphToVectorFilter->SetInput(this->m_Graph);
	graphToVectorFilter->SetOutputDir(this->m_OutputDir);

	//Update filter
	graphToVectorFilter->Update();

	//Update output DS
	this->m_OutputDS = const_cast<OGRDataSourceType*>(graphToVectorFilter->GetOutput());

	//Update layer name
	m_OutputLayerName = cleanedLayerName;

	//Add tile information meta data
	if(m_IsTileProcessing)
	{
		AddMetaData(m_CurrentTile);
	}
	else
	{
		//Create empty tile
		ProcessingTile tile;
		tile.m_Tx = 0;
		tile.m_Ty = 0;

		//Update max tile size
		this->m_MaxTileSizeX = this->m_ImageWidth;
		this->m_MaxTileSizeY = this->m_ImageHeight;

		//add metadata
		AddMetaData(tile);
	}

	//Create filter to simplify
	if(m_IsSimplify)
	{

		auto simplifyVectorFilter = SimplifyFilterType::New();
		simplifyVectorFilter->SetSimplifyFunc(m_SimplifyFunc);
		simplifyVectorFilter->SetInput(graphToVectorFilter->GetOutput());
		simplifyVectorFilter->SetLayerName(otb::obia::cleanedLayerName);
		//FIXME : A modifier pour ne pas casser le pipeline
		//En effet, l'update ne doit se faire que sur le dernier filtre
		//Il faut ajouter un updateOutputInformation
		simplifyVectorFilter->Update();

		//By using New with GDALdataset, the new OGRDatasource get ownership of Dataset.
		//So we have to clone it
		OGRDataSourceType::Pointer tmp = const_cast<OGRDataSourceType*>(simplifyVectorFilter->GetOutput());
		//this->m_OutputDS = OGRDataSourceType::New(VectorOperations::CloneDataset(&tmp->ogr()));
		this->m_OutputDS = const_cast<OGRDataSourceType*>(simplifyVectorFilter->GetOutput());
		//std::cout << "End run filter" << std::endl;
		m_OutputLayerName = reconstructedLayerName;
	}

}
template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::ExtractStabilityMargins()
{
	std::cout << "--------- EXTRACT STABILITY MARGIN --------------" << std::endl;
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    // Compute the number of adjacency layers to extract
    const uint32_t nbAdjacencyLayers = 1;

    // the local serialized margin
    m_MaxNumberOfBytes = 0;

    for(auto& kv : this->m_TileMap)
    {

        uint32_t tx = kv.first % this->m_NumberOfTilesX;
        uint32_t ty = kv.first / this->m_NumberOfTilesX;

		// Load the graph if necessary
		this->ReadGraphIfNecessary(ty, tx, this->m_InputDirectory, this->m_GraphPrefixName);

		std::cout << "Number of nodes = " << this->m_Graph->GetNumberOfNodes() << std::endl;
		// Retrieve the tile
		auto tile = kv.second;

		auto subGraphMap = GraphOperationsType::ExtractStabilityMargin(this->m_Graph,
																		nbAdjacencyLayers,
																		tile,
																		this->m_NumberOfTilesX,
																		this->m_NumberOfTilesY,
																		this->m_ImageWidth
																		  /*this->m_ImageHeight*/);

		std::cout << "Number margin nodes = " << subGraphMap.size() << std::endl;
		m_SerializedStabilityMargin = GraphOperationsType::SerializeStabilityMargin(subGraphMap,
																					this->m_Graph);


		/********************************************************************************************/
		//DEBUG
		//Convert to image
		std::stringstream os;
		os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".tif";
		InputGraphPointerType margeGraph = GraphOperationsType::DeSerializeGraph(m_SerializedStabilityMargin);
		std::cout << "Graph = " << margeGraph->GetNumberOfNodes() << std::endl;
		margeGraph->SetImageWidth(this->m_Graph->GetImageWidth());
		margeGraph->SetImageHeight(this->m_Graph->GetImageHeight());
		margeGraph->SetNumberOfSpectralBands(this->m_Graph->GetNumberOfSpectralBands());
		margeGraph->SetProjectionRef(this->m_Graph->GetProjectionRef());
		ConvertGraphToImage(margeGraph, os.str());
		/********************************************************************************************/


		if(m_MaxNumberOfBytes < m_SerializedStabilityMargin.size())
		{
			m_MaxNumberOfBytes = m_SerializedStabilityMargin.size();
		}

		std::cout << "m_MaxNumberOfBytes = " << m_MaxNumberOfBytes << std::endl;

		if(this->m_TileMap.size() > 1)
		{
			std::stringstream os;
			os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
			GraphOperationsType::WriteSerializedMarginToDisk(m_SerializedStabilityMargin, os.str());
			m_SerializedStabilityMargin.clear();
		}

	}

    mpiConfig->barrier();

    mpiTools->ComputeMax<unsigned long int>(m_MaxNumberOfBytes, MPI_UNSIGNED_LONG);
    std::cout << "m_MaxNumberOfBytes = " << m_MaxNumberOfBytes << std::endl;
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::AggregateStabilityMargins()
{
	std::cout << "--------- AGGREGATE STABILITY MARGIN --------------" << std::endl;
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();


    // Create the shared buffer which will be accessible by other processes.
    uint64_t maxNumberOfElements = this->m_MaxNumberOfTilesPerProcessor * (IntSize + m_MaxNumberOfBytes);
    std::vector< char > sharedBuffer(maxNumberOfElements);

    uint64_t ntile = 0;
    for(auto& kv : this->m_TileMap)
    {

	   uint32_t tx = kv.first % this->m_NumberOfTilesX;
	   uint32_t ty = kv.first / this->m_NumberOfTilesX;

	   uint32_t tile_id = ty*this->m_NumberOfTilesX + tx;

       std::cout << "Share margin of tile " << tx << "/" << ty << std::endl;

		if(this->m_TileMap.size() > 1)
		{
			std::stringstream os;
			os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
			m_SerializedStabilityMargin = GraphOperationsType::ReadSerializedMarginFromDisk(os.str());
		}

		// Move at the right location in the shared buffer.
		uint64_t offset = ntile * (IntSize + m_MaxNumberOfBytes);
		std::cout <<"Offset = " << offset << std::endl;
		std::cout <<"Size   = " << m_SerializedStabilityMargin.size() << std::endl;
		// Write the serialized stablity margin in the shared buffer
		to_stream(sharedBuffer,m_SerializedStabilityMargin,offset);

		// Can release this serialized stability margin
		m_SerializedStabilityMargin.clear();
		m_SerializedStabilityMargin.shrink_to_fit();

		// Increment the number of tiles processed
		ntile++;

    }

    mpiConfig->barrier();

    // Creation of rma window: each processor will have its shared buffer accessible for other processors
    //This is a collectiv call, so we do not need to wait others process to do it
    MPI_Win win;
    MPI_Win_create(&sharedBuffer[0], maxNumberOfElements, CharSize, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    for(auto& kv : this->m_TileMap)
    {

	   uint32_t tx = kv.first % this->m_NumberOfTilesX;
	   uint32_t ty = kv.first / this->m_NumberOfTilesX;
	   std::cout << "Tile " << tx << "/" << ty << " for proc " << mpiConfig->GetMyRank() << std::endl;

	   uint32_t tile_id = ty*this->m_NumberOfTilesX + tx;

		// Retrieve the tile
		auto tile = kv.second;

	   // Creation of the list of serialized stability margins per processor
		std::vector< std::vector<char> > otherSerializedMargins;

		// Retrieve the neighbor tiles
		auto neighborTiles = SpatialTools::EightConnectivity(tile_id, this->m_NumberOfTilesX, this->m_NumberOfTilesY);
		for(unsigned short n = 0; n < 8; n++)
		{
			//std::cout << "Tile id = " << tile_id << ". Neighbor " << n << " is " << neighborTiles[n] << std::endl;
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
				otherSerializedMargins.push_back( std::vector<char>(IntSize + m_MaxNumberOfBytes) );

				// Read rma operation

				//MPI_Win_fence(0, win);
				MPI_Win_lock(MPI_LOCK_SHARED, neighRank, 0, win);

				MPI_Get(&(otherSerializedMargins[otherSerializedMargins.size()-1][0]),
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

		//Read input graph
		this->ReadGraphIfNecessary(ty, tx, this->m_InputDirectory, this->m_GraphPrefixName);

		// Agregate the stability margins to the graph
		for(uint32_t i = 0; i < otherSerializedMargins.size(); i++)
		{
			// Retrieve the serialized margin
			std::vector< char > otherSerializedMargin(0);
			from_stream(otherSerializedMargins[i],otherSerializedMargin);

			otherSerializedMargins[i].clear();
			otherSerializedMargins[i].shrink_to_fit();

			// Deserialize the graph
			auto subGraph = GraphOperationsType::DeSerializeGraph(otherSerializedMargin);

			GraphOperationsType::AggregateGraphs(this->m_Graph, subGraph);

		}

		// Remove duplicated nodes
		auto borderNodeMap = GraphOperationsType::BuildBorderNodesMap(this->m_Graph,
																	  tile,
																	  this->m_NumberOfTilesX,
																	  this->m_NumberOfTilesY,
																	  this->m_ImageWidth);

		GraphOperationsType::RemoveDuplicatedNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth);


		// Update edges
		GraphOperationsType::DetectNewAdjacentNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth, this->m_ImageHeight);

		this->WriteGraphIfNecessary(ty, tx, this->m_TemporaryDirectory, "Graph_With_Marge");

    }

    mpiConfig->barrier();

    // Can release the rma window
    MPI_Win_free(&win);
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::AddMetaData(const ProcessingTile tile)
{
	std::cout << "Adding Metadata to tile " << tile.m_Tx << "_" << tile.m_Ty << std::endl;
	std::ostringstream ss;

	//Compute origin tile X
	ss.clear();
	ss.str("");
	ss << tile.m_Tx*this->m_MaxTileSizeX;
	std::string originX = ss.str();

	//Compute origin tile Y
	ss.str("");
	ss << tile.m_Ty*this->m_MaxTileSizeY;
	std::string originY = ss.str();

	//Add tile size
	ss.str("");
	ss << this->m_MaxTileSizeX;
	std::string tileSizeX = ss.str();

	//Add tile size
	ss.str("");
	ss << this->m_MaxTileSizeY;
	std::string tileSizeY = ss.str();

	//Add origin X image
	ss.str("");
	ss << this->m_Graph->GetOriginX();
	std::string imageOriginX = ss.str();

	//Add origin Y image
	ss.str("");
	ss << this->m_Graph->GetOriginY();
	std::string imageOriginY = ss.str();

	this->m_OutputDS->ogr().SetMetadataItem("OriginTileX"	, originX.c_str());
	this->m_OutputDS->ogr().SetMetadataItem("OriginTileY"	, originY.c_str());
	this->m_OutputDS->ogr().SetMetadataItem("TileSizeX"  	, tileSizeX.c_str());
	this->m_OutputDS->ogr().SetMetadataItem("TileSizeY" 	, tileSizeY.c_str());
	this->m_OutputDS->ogr().SetMetadataItem("ImageOriginX"  , imageOriginX.c_str());
	this->m_OutputDS->ogr().SetMetadataItem("ImageOriginY"  , imageOriginY.c_str());
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::RemovePolygonsOutsideTile(const ProcessingTile& tile)
{
	//Create a polygon equivalent to the tile
	OGRPolygon* tilePolygon = CreateTilePolygon(tile);

	//Loop accross features
	OGRLayerType layer = this->m_OutputDS->GetLayer(m_OutputLayerName);
	unsigned int nbFeatures = layer.GetFeatureCount(true);

	for(int fId = 0; fId < nbFeatures; ++fId)
	{
		//Current feature
		OGRFeatureType curFeature = layer.GetFeature(fId);

		//Get geometry
		const OGRGeometry* curGeom = curFeature.GetGeometry();

		if(!curGeom->Intersects(tilePolygon))
		{
			//If the geometry do not intersects polygon tile, so this a feature belonging to stability margin
			//Destroy the feature
			layer.DeleteFeature(fId);
		}
	}

	delete tilePolygon;
}

template< class TInputGraph, class TSimplifyFunc >
OGRPolygon*
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::CreateTilePolygon(const ProcessingTile& tile)
{
	int xStart = tile.m_Frame.GetIndex().GetElement(0);
	int yStart = tile.m_Frame.GetIndex().GetElement(1);
	int xSize = tile.m_Frame.GetSize(0);
	int ySize = tile.m_Frame.GetSize(1);

	OGRPolygon* tilePolygon = new OGRPolygon();
	OGRLinearRing* extRing = new OGRLinearRing();

	//Add all points
	extRing->addPoint(xStart, yStart);
	extRing->addPoint(xStart + xSize - 1, yStart);
	extRing->addPoint(xStart + xSize - 1, yStart + ySize - 1);
	extRing->addPoint(xStart, yStart + ySize - 1);
	extRing->addPoint(xStart, yStart);

	//Add ring
	tilePolygon->addRingDirectly(extRing);

	return tilePolygon;
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::WriteFeatures(const ProcessingTile& tile)
{
	std::cout << "WRITING FEATURES FOR TILE " << tile.m_Ty << "_" << tile.m_Tx << std::endl;
	//Create filename
	std::stringstream ss;
	ss << this->m_OutputDir << "Reconstructed_Polygons_" << tile.m_Tx << "_" << tile.m_Ty << ".gml";

	VectorOperations::WriteOGRDataSource(this->m_OutputDS, ss.str(), "");

//	//Get output
//	GDALDriver*  poDriverGml = VectorOperations::InitializeGDALDriver("GML");
//	GDALDataset* polyGDs = poDriverGml->Create(ss.str().c_str(), 0, 0, 0, GDT_Unknown, NULL );
//
//	std::cout << "Number of layer = " << this->m_OutputDS->GetLayersCount() << std::endl;
//	//Write all layers
//	for(unsigned int layerId = 0; layerId < this->m_OutputDS->ogr().GetLayerCount(); ++layerId)
//	{
//		OGRLayerType curLayer =  this->m_OutputDS->GetLayer(layerId);
//		std::cout << "Number of polygon to write = " << curLayer.GetFeatureCount(true) << std::endl;
//		//std::cout << "Layer = " << curLayer.GetName() << std::endl;
//		polyGDs->CopyLayer(&(curLayer.ogr()), curLayer.GetName().c_str());
//	}
//
//
//	//Clean memory
//	GDALClose(polyGDs);
}

template< class TInputGraph, class TSimplifyFunc >
void
LSPolygonizeScheduler<TInputGraph, TSimplifyFunc>
::ConvertGraphToImage(InputGraphPointerType inputGraph, std::string filename)
{
//	using LabelPixelType = unsigned int;
//	using LabelImageType = otb::Image< LabelPixelType, 2 >;
//	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<TInputGraph, LabelImageType>;
//	using RGBPixelType = itk::RGBPixel<unsigned char>;
//	using RGBImageType = otb::Image<RGBPixelType, 2>;
//	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
//	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
//	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;
//
//
//	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
//	auto labelToRGBFilter = LabelToRGBFilterType::New();
//	auto rgbWriter = RGBWriterType::New();
//
//	rgbWriter->SetFileName(filename);
//	graphToLabelFilter->SetInput(inputGraph);
//	auto fillHoleFilter = FillholeFilterType::New();
//	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
//	labelToRGBFilter->SetInput(fillHoleFilter->GetOutput());
//	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
//	rgbWriter->Update();


	using LabelPixelType = unsigned int;
	using LabelImageType = otb::Image< LabelPixelType, 2 >;
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<TInputGraph, LabelImageType>;
	using WriterType = otb::ImageFileWriter< LabelImageType>;

	/** FOR COLOR IMAGE
	using RGBPixelType = itk::RGBPixel<unsigned char>;
	using RGBImageType = otb::Image<RGBPixelType, 2>;
	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
	*/
	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;

	//Output name
	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	auto grayWriter = WriterType::New();
	auto fillHoleFilter = FillholeFilterType::New();

	grayWriter->SetFileName(filename);
	graphToLabelFilter->SetInput(inputGraph);
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	grayWriter->SetInput(fillHoleFilter->GetOutput());
	grayWriter->Update();
}
} // end of namespace obia
} // end of namespace otb

#endif
