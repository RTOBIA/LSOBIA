#ifndef otbObiaLSLSSmallRegionsMergingFilter_txx
#define otbObiaLSLSSmallRegionsMergingFilter_txx

#include "otbObiaGraphOperations.h"
#include "otbObiaSmallRegionsMergingFilter.h"
#include "otbObiaLSSmallRegionsMergingScheduler.h"
#include "otbObiaMPITools.h"
#include "otbImage.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbObiaGraphToLabelImageFilter.h"
#include "otbImageFileWriter.h"
#include <unordered_set>

namespace otb
{
namespace obia
{

template< typename TGraph >
LSSmallRegionsMergingScheduler<TGraph>
::LSSmallRegionsMergingScheduler() :
m_MinimalSurface(1),
m_MaxNumberOfBytes(0)
{

}

template< typename TGraph >
LSSmallRegionsMergingScheduler<TGraph>
::~LSSmallRegionsMergingScheduler()
{
}


template< class TGraph >
itk::SmartPointer<typename LSSmallRegionsMergingScheduler<TGraph>::SRMFilterType> LSSmallRegionsMergingScheduler<TGraph>
::CreateFilter()
{
	// Merging filter
	auto srmFilter = SRMFilterType::New();

	//Create the 3 classes
	auto mergingFunc = new SRMMergingCost<float, GraphType>();
	mergingFunc->SetMinimalSurface(this->m_MinimalSurface);

	auto heuristicFunc = new SRMHeuristic<GraphType>();//::New();
	heuristicFunc->SetGraph(this->m_Graph);
	heuristicFunc->SetMinimalSurface(this->m_MinimalSurface);

	auto updateFunc = new SRMUpdateAttribute<GraphType>();

	srmFilter->SetMergingCostFunc(mergingFunc);
	srmFilter->SetHeuristicFunc(heuristicFunc);
	srmFilter->SetUpdateAttributeFunc(updateFunc);
	srmFilter->SetMaxNumberOfIterations(this->m_NumberOfIterations);
	srmFilter->SetInput(this->m_Graph);
	return srmFilter;

}

//
//
//template< class TGraph >
//itk::SmartPointer<typename LSSmallRegionsMergingScheduler<TGraph>::SRMFilterType> LSSmallRegionsMergingScheduler<TGraph>
//::CreateFilter()
//{
//	// Merging filter
//	auto srmFilter = SRMFilterType::New();
//
//	srmFilter->SetMinimalSurface(this->m_MinimalSurface);
//	srmFilter->SetNumberOfIterations(this->m_NumberOfIterations);
//	srmFilter->SetInput(this->m_Graph);
//	return srmFilter;
//
//}

template< typename TGraph >
void LSSmallRegionsMergingScheduler<TGraph>
::AggregateGraph()
{
	std::cout << "AGGREGATE FINAL GRAPH" << std::endl;
	//TODO
	std::stringstream os;
	int tx;
	int ty;

}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::GenerateData()
{
    std::cout << "########### Nombre tuile " << this->m_TileMap.size() << std::endl;
    auto mpiConfig = otb::MPIConfig::Instance();
    mpiConfig->barrier();

    if(this->m_TilesPerProcessor.size() == 1 && this->m_TilesPerProcessor.begin()->second.size() == 1)
    {
        std::cout << "No Tiling Execution" << std::endl;
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
        std::cout << "Tiling Execution : TileMap Size : "  << this->m_TileMap.size() << " et  Nombre proc : " << this->m_TilesPerProcessor.size()  << std::endl;
        TilingExecution();
    }
}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::ComputePaddingValue()
{

}

template< typename TGraph >
uint32_t
LSSmallRegionsMergingScheduler<TGraph>
::ComputeMaxDepth()
{
    /**The idea to compute the padding value is the following one:
     *  Explore the graph to compute the max required adjency layer in order to merge all smalls regions
     *  So, for all nodes in the tile's border, we look for nodes whith small area, and then look for adjency node*/

    //Max depth
    uint32_t maxDepth = std::numeric_limits<uint32_t>::min();

    //For each tile
    uint32_t tid = 0;
    for(auto& kv : this->m_TileMap)
    {
        uint32_t tx = kv.first % this->m_NumberOfTilesX;
        uint32_t ty = kv.first / this->m_NumberOfTilesX;

        //std::cout << "Compute Max Depth for " << ty << "_" << tx << std::endl;

        // Retrieve the tile by reference since it will be modified.
        auto& tile = kv.second;

        //Get the graph associated to this tile
        this->ReadGraphIfNecessary(ty, tx);
        /*std::stringstream os;
        os << this->m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
        auto graph = GraphOperationsType::ReadGraphFromDisk(os.str());
        graph->SetImageWidth(this->m_ImageWidth);
        graph->SetImageHeight(this->m_ImageHeight);
        graph->SetNumberOfSpectralBands(this->m_NumberOfSpectralBands);*/

        /**Get Borders Node*/
        // Retrieve the list of the nodes located at the borders of the tile
        auto borderNodes = GraphOperationsType::GetListOfBorderNodes(this->m_Graph, tile, this->m_NumberOfTilesX,
                                                                     this->m_NumberOfTilesY, this->m_ImageWidth);

        /** Extract small regions nodes*/
        std::vector<NodeType*> smallRegionsNodes = ExtractSmallRegionsNodes(borderNodes);
        /*std::cout << "Image Width " << this->m_ImageWidth << " / nTx = " << this->m_NumberOfTilesX << " //" << this->m_NumberOfTilesY << std::endl;
        std::cout << "Nombre noeuds graph " << graph->GetNumberOfNodes() << " pour la tuile : " << ty << "_" << tx <<  std::endl;
        std::cout << "Nombre noeuds en bordures: "<< borderNodes.size() <<  " pour la tuile : " << ty << "_" << tx <<  std::endl;*/
         //std::cout << "Nombre de Small Nodes :  " << smallRegionsNodes.size() << " pour la tuile : " << ty << "_" << tx <<  std::endl;

        /** Compute depth for these node and get the max*/
        std::unordered_map< NodeType*, uint32_t > extractedNodes;
        for(auto& n : smallRegionsNodes)
        {
            RecursiveSmallRegionsDepthBreadthFirstExploration(n, extractedNodes, this->m_Graph, 0, n->m_Attributes.m_Area);
        }

        /**Get the maximum in the map*/
        for(auto& n : extractedNodes)
        {

            if(n.second >= maxDepth)
            {
                maxDepth = n.second;
            }
        }
        //std::cout << "Depth for " << ty << "_" << tx << " = " << maxDepth << std::endl;
    }

    //std::cout <<  "Max depth: " << maxDepth << std::endl;
    /**TODO : maybe use MPI to ensure that the padding value will be the same for all tiles and processor*/
    //this->m_PaddingValue = 0;
    return maxDepth;
}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::NoTilingExecution()
{
    /**No tiling, so just run the SRM Filter on the graph*/
    //Execute the small region filter
    bool merge_over = false;
    std::cout << "Surface : "<< this->m_MinimalSurface << std::endl;
    int i = 0;
    while(!merge_over)
    {
        auto srmFilter = CreateFilter();
        srmFilter->Update();
        this->m_Graph = srmFilter->GetOutput();
        //this->m_Graph->GraftGraphByMove(srmFilter->GetOutput());
        merge_over = srmFilter->GetMergingOver();
        //merge_over = SRMFilter->GetMergeOver();
        std::cout << "Nombre Noeuds after: " << srmFilter->GetInput()->GetNumberOfNodes() << std::endl;
        i++;
        //if(i > 3)
        //    merge_over = true;
    }

    //Once all small regions merged, get the graph
    std::cout << "Nombre Noeuds : " << this->m_Graph->GetNumberOfNodes() << std::endl;
}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::TilingExecution()
{
	/**If we have one tile per processor, then the graph is already initialized (this->m_Graph) from the previous filter
	 * We have to keep this graph and work with it (update, merging, etc ...)*/
	auto mpiConfig = MPIConfig::Instance();
	int localFusion = 0;

	bool merge_over = false;
	uint32_t cur_it = 0;

	int maxNumberOfIterations = 100;

	while(!merge_over)
	{

		int globalFusion = 0;
		std::cout << "Tiling Execution " << cur_it << " for " << MPIConfig::Instance()->GetMyRank() << std::endl;

		/**Compute stability margins*/
		std::cout << "------- EXTRACT -------------" << std::endl;
		ExtractStabilityMargins();

		/** Aggregate*/
		std::cout << "------- AGGREGATE ------------- " << std::endl;
		AggregateStabilityMargins();

		unsigned long int accumulatedMemory = this->m_AvailableMemory + 1;

		/** Run merging small regions*/
		localFusion = PartialFusion(accumulatedMemory, globalFusion);
		/*if(this->m_TileMap.size() == 1){
			std::stringstream os;
			os << this->m_TemporaryDirectory << "AFTER_" << mpiConfig->GetMyRank() << ".dat";
			GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());
		}*/

		if(localFusion == 0)
		{
			std::cout <<"Local fusion = 0" << std::endl;
			merge_over = true;
		}

		if(globalFusion != 0){
			merge_over = false;
		}

		if(cur_it > maxNumberOfIterations){
			merge_over = true;
		}

		cur_it++;

	}

	/*if(m_AggregateGraph)
	{
		AggregateGraph();
	}*/
}


template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::ExtractStabilityMargins()
{
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    // Compute the number of adjacency layers to extract
   // const uint32_t nbAdjacencyLayers = (m_NumberOfIterations < 0) ? ComputeMaxDepth() : m_NumberOfIterations;//ComputeMaxDepth();
   // const uint32_t nbAdjacencyLayers = pow(2, m_NumberOfIterations + 1) - 2;
    const uint32_t nbAdjacencyLayers = pow(2, 1 + 1) - 2;
    std::cout << "Number of adjacency layers : " << nbAdjacencyLayers << std::endl;
    // the local serialized margin
    std::vector< char > serializedMargin;
    m_MaxNumberOfBytes = 0;

    for(auto& kv : this->m_TileMap)
    {
        uint32_t tx = kv.first % this->m_NumberOfTilesX;
        uint32_t ty = kv.first / this->m_NumberOfTilesX;

        // Retrieve the tile by reference since it will be modified.
        auto& tile = kv.second;

        // Load the graph if necessary
        this->ReadGraphIfNecessary(ty, tx);

        //std::cout << "Graph "  << ty << "_" << tx  <<  " read : " << this->m_Graph->GetNumberOfNodes() << std::endl;
		std::cout << "------ TILE : " << ty << "_" << tx << " for processor : " << mpiConfig->GetMyRank()<< std::endl;

        auto subGraphMap = GraphOperationsType::ExtractStabilityMargin(this->m_Graph,
                                                                       nbAdjacencyLayers,
                                                                       tile,
                                                                       this->m_NumberOfTilesX,
                                                                       this->m_NumberOfTilesY,
                                                                       this->m_ImageWidth
                                                                       /*this->m_ImageHeight*/);

        std::cout << "Margin Graph (number of nodes ) : " << subGraphMap.size() << std::endl;

        m_SerializedStabilityMargin = GraphOperationsType::SerializeStabilityMargin(subGraphMap,
                                                                                    this->m_Graph);

        if(m_MaxNumberOfBytes < m_SerializedStabilityMargin.size())
        {
            m_MaxNumberOfBytes = m_SerializedStabilityMargin.size();
        }

        //If we have several tiles per processor, we have to store the graph on the disk
        if(this->m_TileMap.size() > 1)
        {
            std::stringstream os;
            os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
            GraphOperationsType::WriteSerializedMarginToDisk(m_SerializedStabilityMargin, os.str());
            m_SerializedStabilityMargin.clear();
        }
    }
    std::cout << "WAIT ALL..Extracting current proc is " << mpiConfig->GetMyRank() << std::endl;
    mpiConfig->barrier();

    mpiTools->ComputeMax<unsigned long int>(m_MaxNumberOfBytes, MPI_UNSIGNED_LONG);
    std::cout << "Max Bytes = " << m_MaxNumberOfBytes << std::endl;
}

template< typename TGraph >
std::vector< char >
LSSmallRegionsMergingScheduler<TGraph>
::ShareStabilityMargins()
{
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

        std::cout << "------ SHARE TILE : " << ty << "_" << tx << " for processor : " << mpiConfig->GetMyRank()<< std::endl;
        // Retrieve the tile by reference since it will be modified.
        auto& tile = kv.second;

        if(this->m_TileMap.size() > 1)
        {
            std::stringstream os;
            os << this->m_TemporaryDirectory << "MarginGraph_" << ty << "_" << tx << ".dat";
            m_SerializedStabilityMargin = GraphOperationsType::ReadSerializedMarginFromDisk(os.str());
        }

        // Move at the right location in the shared buffer.
        uint64_t offset = ntile * (IntSize + m_MaxNumberOfBytes);

        // Write the serialized stablity margin in the shared buffer
	to_stream(sharedBuffer,m_SerializedStabilityMargin,offset);

        // Can release this serialized stability margin
        m_SerializedStabilityMargin.clear();
        m_SerializedStabilityMargin.shrink_to_fit();

        // Increment the number of tiles processed
        ntile++;

    } // end for(uint32_t ty = 0; ty < this->m_NumberOfTilesY; ty++)

    //Wait ALL
    mpiConfig->barrier();

    std::cout << "WAIT ALL...Sharing..., current proc is " << mpiConfig->GetMyRank() << std::endl;
    return sharedBuffer;
}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::AggregateStabilityMargins()
{
	auto mpiConfig = MPIConfig::Instance();
	auto mpiTools = MPITools::Instance();

	// Create the shared buffer which will be accessible by other processes.
	std::vector< char > sharedBuffer = ShareStabilityMargins();

	MPI_Win win;
	MPI_Win_create(&sharedBuffer[0], sharedBuffer.size(), CharSize, MPI_INFO_NULL, MPI_COMM_WORLD, &win);


	uint64_t ntile = 0;
	for(auto& kv : this->m_TileMap)
	{
		uint32_t tx = kv.first % this->m_NumberOfTilesX;
		uint32_t ty = kv.first / this->m_NumberOfTilesX;

		std::cout << "------ AGGREGATE TILE : " << ty << "_" << tx << " for processor : " << mpiConfig->GetMyRank()<< std::endl;

		// Retrieve the tile by reference since it will be modified.
		auto& tile = kv.second;
		// Creation of the list of serialized stability margins per processor
		std::vector< std::vector<char> > otherSerializedMargins;

		// Retrieve the neighbor tiles
		uint32_t tile_id = ty*this->m_NumberOfTilesX + tx;
		auto neighborTiles = SpatialTools::EightConnectivity(tile_id, this->m_NumberOfTilesX, this->m_NumberOfTilesY);
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

		//Read the current graph
		this->ReadGraphIfNecessary(ty, tx);

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

			std::cout << "Sub Graph : " << subGraph->GetNumberOfNodes() << " for tile "   << ty << "_" << tx << std::endl;

			std::cout << "Adding " << this->m_Graph->GetNumberOfNodes()
					  << " and " << subGraph->GetNumberOfNodes() << " nodes " << std::endl;
			//Add the sub graph to the graph
			GraphOperationsType::AggregateGraphs(this->m_Graph, subGraph);

			//Free memory
			subGraph->Reset();
		}

		// Remove duplicated nodes
		auto borderNodeMap = GraphOperationsType::BuildBorderNodesMap(this->m_Graph,
																	  tile,
																	  this->m_NumberOfTilesX,
																	  this->m_NumberOfTilesY,
																	  this->m_ImageWidth);

		uint32_t nb_nodes_before = this->m_Graph->GetNumberOfNodes();
		GraphOperationsType::RemoveDuplicatedNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth);

		// Update edges
		GraphOperationsType::DetectNewAdjacentNodes(borderNodeMap, this->m_Graph, this->m_ImageWidth, this->m_ImageHeight);

		//Clear borderNodeMap
		borderNodeMap.clear();

		this->WriteGraphIfNecessary(ty, tx);
		ntile++;
	}

	std::cout << "WAIT ALL AGGREGATE, current proc is " << mpiConfig->GetMyRank() << std::endl;

	//Wait All
	mpiConfig->barrier();

	// Can release the rma window
	MPI_Win_free(&win);
}


template< typename TGraph >
int
LSSmallRegionsMergingScheduler<TGraph>
::PartialFusion(unsigned long int& accumulatedMemory,  int& globalFusion)
{
    std::cout << "--------- FUSION (minimal surface : " << m_MinimalSurface << ")-----------" << std::endl;

    /**Parcours des tuiles*/
    auto mpiConfig = MPIConfig::Instance();
    auto mpiTools = MPITools::Instance();

    accumulatedMemory = 0;
    int localFusion  = 0;

    for(auto& kv : this->m_TileMap)
    {
        uint32_t tx = kv.first % this->m_NumberOfTilesX;
        uint32_t ty = kv.first / this->m_NumberOfTilesX;

        // Retrieve the tile by reference since it will be modified.
        auto& tile = kv.second;
        std::cout << "---- FUSION TUILE " << ty << "_" << tx << " for processor " << mpiConfig->GetMyRank() <<  std::endl;
        fflush(stdout);

        //Get the graph associated to this tile
        this->ReadGraphIfNecessary(ty, tx);
        std::cout << "Graph : image Width = " << this->m_Graph->GetImageWidth() << "/ Nombre noeuds = " << this->m_Graph->GetNumberOfNodes() <<std::endl;


        //Execute the small region filter
        auto srmFilter  = CreateFilter();
        srmFilter->Update();

        std::cout << "After update" << std::endl;
        //srmFilter->GetOutput()->CopyGraph(&(*this->m_Graph));
        this->m_Graph = srmFilter->GetOutput();

        //Remove unstable
        std::cout <<"Remove Unstable" << std::endl;
        GraphOperationsType::RemoveUnstableNodes(this->m_Graph,
                                                 tile,
                                                 this->m_ImageWidth);


        std::cout << "Nombre after remove : " << this->m_Graph->GetNumberOfNodes() << std::endl;
        localFusion += (srmFilter->GetMergingOver()) ? 0 : 1;
        this->WriteGraphIfNecessary(ty, tx);

        std::cout << "----  END FUSION TUILE " << ty << "_" << tx << " for processor " << mpiConfig->GetMyRank() <<  std::endl;
    }

    //Update globalFusion
    globalFusion = localFusion;


    //std::cout << "End of processor " << mpiConfig->GetMyRank() << std::endl;
    fflush(stdout);
    //Wait all procs finished
    std::cout << "WAIT AFTER PARTIEL : " << mpiConfig->GetMyRank() << std::endl;
    mpiConfig->barrier();

    /* Compute the accumulated memory */
    mpiTools->Accumulate(accumulatedMemory,  MPI_UNSIGNED_LONG);

    /* Compute global fusion*/
    mpiTools->ComputeMax(globalFusion, MPI_INTEGER);

    return localFusion;

}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::ReconditionGraph()
 {
    auto mpiConfig = MPIConfig::Instance();
    //Number of nodes
    uint32_t nbNodes = this->m_Graph->GetNumberOfNodes();

    //Looping across nodes
    //for(auto nodeIt = outputGraph->Begin(); nodeIt != outputGraph->End(); nodeIt++)
    for(uint32_t k = 0; k < nbNodes; ++k)
    {
        //std::cout << "Noeud " << k << std::endl;
        //NodeType* nodeIt = this->m_Graph->GetNodeAt(k);
        //nodeIt->m_Valid = true;
        /*if(nodeIt->m_Attributes.m_Area < m_MinimalSurface){
            std::cout << "Graph pas complet " << std::endl;
            exit(1);
        }*/
    }

    //Write graph
    /*std::stringstream os;
    os << this->m_TemporaryDirectory << "Graph_" << mpiConfig->GetMyRank()<< ".dat";
    GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());

    //Read graph
    this->m_Graph = GraphOperationsType::ReadGraphFromDisk(os.str());*/

 }


template< typename TGraph >
std::vector< typename TGraph::NodeType* >
LSSmallRegionsMergingScheduler<TGraph>
::ExtractSmallRegionsNodes(std::vector< typename TGraph::NodeType* > borderNodes)
{
    std::vector< NodeType* > smallRegionsNodes;

    // For each node we visit its adjacency layers using the DFS graph algorithm.
    for(auto& n : borderNodes)
    {
        //If node's area is smaller than the minimal surface, then add to the list
        if(n->m_Attributes.m_Area <= m_MinimalSurface)
        {
            smallRegionsNodes.push_back(n);
        }
    }

    return smallRegionsNodes;
}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::RecursiveSmallRegionsDepthBreadthFirstExploration(typename TGraph::NodeType* node,
                                                                       std::unordered_map< typename TGraph::NodeType*, uint32_t >& visitedNodes,
                                                    typename TGraph::Pointer graph,
                                                    const uint32_t currentNumberOfAdjacencyLayers,
                                                    const uint32_t currentArea)
{
    //std::cout << "Explore node " << node->m_Id << " with area " << currentArea <<  " and depth "  << currentNumberOfAdjacencyLayers << std::endl;


    //Vérifier aussi que les noeuds explorés n'appartiennent pas à la bordure...
    // Try to retrieve the depth exploration of this node
    auto isHere = visitedNodes.find(node);
    // The node has been already visited
    if(isHere != visitedNodes.end())
    {
        //std::cout << " Explore " << node->m_Id << " again ? " << std::endl;
        // If the current depth exploration is lower than the node 's depth exploration
        // then we update its depth with this one and we re-explore this node.
        if(currentNumberOfAdjacencyLayers < isHere->second)
        {
            // Update the depth exploration
            isHere->second = currentNumberOfAdjacencyLayers;

            // Loop over its edges to explore recursively its adjacent nodes.
            for(auto edg : node->m_Edges)
            {
                //If the node has a small area, then explore it

                //If the cumulative area is smaller than the minimal surface, then we keep exploring
                uint32_t cumulativeArea = currentArea + graph->GetNodeAt(edg.m_TargetId)->m_Attributes.m_Area;

                if(cumulativeArea <= m_MinimalSurface)
                {
                    RecursiveSmallRegionsDepthBreadthFirstExploration(graph->GetNodeAt(edg.m_TargetId),
                                                                      visitedNodes, graph, currentNumberOfAdjacencyLayers + 1,
                                                                      cumulativeArea);
                }

            } // end for for(auto edg : node->m_Edges)

        } // end if(currentNumberOfAdjacencyLayers < isHere->second)
        else
        {
            //std::cout << "  Damn NO!" << std::endl;
            // Second stopping criterion: this node has been fully explored.
            return;
        }
    } // end if
    else
    {
        //std::cout << "First time for " << node->m_Id << std::endl;
        // This node is explored for the first time.
        visitedNodes[node] = currentNumberOfAdjacencyLayers;

        // Loop over its edges to explore recursively its adjacent nodes.
        for(auto edg : node->m_Edges)
        {

            //If the cumulative area is smaller than the minimal surface, then we keep exploring
            uint32_t cumulativeArea = currentArea + graph->GetNodeAt(edg.m_TargetId)->m_Attributes.m_Area;
            //std::cout << "cumulativeArea for " << node->m_Id <<  " and node " << graph->GetNodeAt(edg.m_TargetId)->m_Id << " = " << cumulativeArea << std::endl;
            //std::cout << "De we Keep exploring ?(" << cumulativeArea << " <= " << m_MinimalSurface << ")" << std::endl;
            if(cumulativeArea <= m_MinimalSurface)
            {

                //Update adjency layers, and explore this new node with the current area
                RecursiveSmallRegionsDepthBreadthFirstExploration(graph->GetNodeAt(edg.m_TargetId),
                                                                  visitedNodes, graph, currentNumberOfAdjacencyLayers + 1,
                                                                  cumulativeArea);
            }

        } // end for(auto edg : node->m_Edges)

    } // end else

} // end else

template< typename TGraph >
bool
LSSmallRegionsMergingScheduler<TGraph>
::IsBorderNode(NodeType* node, const std::vector< typename GraphOperations<TGraph>::NodeType* > borderNodes)
 {
    //std::cout << "Is Node " << node->GetFirstPixelCoords() << " is border ? " << std::endl;
    for(auto& n : borderNodes)
    {
        //std::cout <<" Border : " << n->GetFirstPixelCoords()  << std::endl;
        if(n->GetFirstPixelCoords() == node->GetFirstPixelCoords())
        {
            return true;
        }
    }

    return false;
 }

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::CreateOutput()
{
    if(this->m_TilesPerProcessor.size() == 1 && this->m_TilesPerProcessor.begin()->second.size() == 1)
    {
        std::stringstream os;
        os << this->m_TemporaryDirectory << "Graph_SRM.dat";
        GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());
    }
    for(auto& kv : this->m_TileMap)
    {
            uint32_t tx = kv.first % this->m_NumberOfTilesX;
            uint32_t ty = kv.first / this->m_NumberOfTilesX;

            std::stringstream os;
            os << this->m_TemporaryDirectory << "Graph_SRM_" << ty << "_" << tx << ".dat";
            GraphOperationsType::WriteGraphToDisk(this->m_Graph, os.str());

            //this->ReadGraphIfNecessary(ty, tx);
            if(m_WriteImage)
            {
            	ConvertToImage(ty, tx);
            }
    }

}

template< typename TGraph >
void
LSSmallRegionsMergingScheduler<TGraph>
::ConvertToImage(uint32_t ty, uint32_t tx)
{

	std::cout << "WRITE LABEL IMAGE FOR " << ty <<"_"<< tx << std::endl;
	if(this->m_Graph == nullptr)
	{
		std::cout << "Graph " << ty << "_" << tx << " is not set..." << std::endl;
		exit(EXIT_FAILURE);
	}
	using LabelPixelType = unsigned int;
	using LabelImageType = otb::Image< LabelPixelType, 2 >;
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<GraphType, LabelImageType>;
	using WriterType = otb::ImageFileWriter< LabelImageType>;

	//FOR COLOR IMAGE
	//      	using RGBPixelType = itk::RGBPixel<unsigned char>;
	//      	using RGBImageType = otb::Image<RGBPixelType, 2>;
	//      	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
	//      	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;

	using FillholeFilterType           = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;

	//Output name
	std::stringstream os;
	os << this->m_OutputDir << "Graph_SRM_" << ty << "_" << tx << ".tif";

	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	auto grayWriter = WriterType::New();
	auto fillHoleFilter = FillholeFilterType::New();
	//auto labelRGB = LabelToRGBFilterType::New();
	//auto rgbWriter = RGBWriterType::New();

	graphToLabelFilter->SetInput(this->m_Graph);
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	//labelRGB->SetInput(fillHoleFilter->GetOutput());
	//rgbWriter->SetInput(labelRGB->GetOutput());
	//rgbWriter->SetFileName(os.str());
	//rgbWriter->Update();
	//
	grayWriter->SetFileName(os.str());
	grayWriter->SetInput(fillHoleFilter->GetOutput());
	grayWriter->Update();
}




}//End obia
}//End otb

#endif
