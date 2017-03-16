#include <iostream>

#include "otbObiaBaatzSegmentationFilter.txx"
#include "otbObiaBaatzGraphToSRMGraphFilter.txx"
#include "otbObiaImageToBaatzGraphFilter.txx"
#include "otbObiaLSMeanShiftScheduler.txx"
#include "otbObiaSmallRegionsMergingFilter.txx"



int otbObiaSplitGraph(int argc, char * argv[])
{
	if(argc < 4)
	{

	  std::cerr << "Usage " << argv[0] << ": \n"
			<< "Argument 1: " << "[path to the input graph]\n"
			<< "Argument 2: " << "[temporary directory to store intermediate files for a node]\n"
			<< "Argument 3: " << "[Tile Size X]\n"
			<< "Argument 4: " << "[Tile Size Y]"
			<< std::endl;

	  return 1;
	}

	using NodeType                     = otb::obia::Node< otb::obia::BaatzNodeAttribute, otb::obia::BaatzEdgeAttribute >;
	using EdgeType                     = typename NodeType::EdgeType;
	using InputGraphType              = otb::obia::Graph< NodeType >;
	using InputGraphPointerType       = typename InputGraphType::Pointer;
	using GraphOperationsType         = otb::obia::GraphOperations<InputGraphType>;

	// Input parameters
	const std::string filename = argv[1];
	const std::string tmpDir = argv[2];
	const uint32_t tileSizeX = atoi(argv[3]);
	const uint32_t tileSizeY = atoi(argv[4]);


	using InputImageType              = otb::VectorImage<float, 2>;
	using LabelPixelType              = unsigned int;
	using SmallRegionMergingFilter	= otb::obia::SmallRegionsMergingFilter<InputGraphType>;

	//Lecture gu graphe
	std::cout << "Tile Size X " << tileSizeX << "/ Tile Size Y " << tileSizeY << std::endl;

	std::stringstream os;
	os << filename;

	auto graph = GraphOperationsType::ReadGraphFromDisk(os.str());

	//Compute number of tiles X and Y
	std::vector<otb::obia::ProcessingTile> tilesList;
	
	uint32_t imageWidth  = graph->GetImageWidth();
	uint32_t imageHeight = graph->GetImageHeight();

	std::cout << "Image width " << imageWidth << std::endl;
	uint32_t numberOfTilesX = ceil(imageWidth / tileSizeX);
	uint32_t numberOfTilesY = ceil(imageHeight/ tileSizeY);
	uint32_t tid = 0;

	for(uint32_t ty = 0; ty < numberOfTilesY; ty++)
	{
		for(uint32_t tx = 0; tx < numberOfTilesX; tx++)
		{

			std::cout << "Creating tile " << ty << "_" << tx << std::endl;

			otb::obia::ProcessingTile tile;

		  	tile.m_Tx = tx;
		  	tile.m_Ty = ty;
		  	tile.m_Frame.SetIndex(0, tx * tileSizeX);
		  	tile.m_Frame.SetIndex(1, ty * tileSizeY);
		  	tile.m_Frame.SetSize(0, std::min((long int)tileSizeX, (long int)imageWidth - tile.m_Frame.GetIndex(0)));
		  	tile.m_Frame.SetSize(1, std::min((long int)tileSizeY, (long int)imageHeight - tile.m_Frame.GetIndex(1)));
			tile.m_MarginValues.fill(0);


		      tilesList.push_back(tile);

		} // end if (mpiTools->IsMyTurn(tid))

		    tid++;

	} // end for(uint32_t tx = 0; tx < nbTilesX; tx++)

	//Release the graph
	graph->Reset();

	//For each tile, remove nodes outside this tile
	for(uint32_t k = 0; k < tilesList.size(); ++k)
	{
		otb::obia::ProcessingTile tile = tilesList[k];

		std::cout << "Tile " << tile.m_Ty << "_" << tile.m_Tx << std::endl;
		std::cout << "Size margin 0 = " << tile.m_MarginValues[0] << std::endl;
		std::cout << "Size margin 1 = " << tile.m_MarginValues[0] << std::endl;
		std::cout << "Size margin 2 = " << tile.m_MarginValues[0] << std::endl;
		std::cout << "Size margin 3 = " << tile.m_MarginValues[0] << std::endl;
		std::cout << "Read graph from disk : " << os.str() << std::endl;
		//Create a new graph
		auto subGraph =  GraphOperationsType::ReadGraphFromDisk(os.str());

		std::cout << "Nombre noeud : " << subGraph->GetNumberOfNodes() << std::endl;
		std::cout << "Noeud 1 " << subGraph->GetNodeAt(0)->m_Attributes.m_Area << std::endl;
		std::cout << "Image width " << subGraph->GetImageWidth() << std::endl;

		//remove nodes outside the tile
		GraphOperationsType::RemoveUnstableNodes(subGraph, tile, subGraph->GetImageWidth());

		//Write Graph to disk
		std::stringstream oso;
		oso << tmpDir <<"/"<< "Splitted_Graph_" << tile.m_Ty << "_" << tile.m_Tx << ".dat" ;
		GraphOperationsType::WriteGraphToDisk(subGraph, oso.str());

		subGraph->Reset();
	}

	std::cout << "Sucess" << std::endl;

	  return EXIT_SUCCESS;
}
