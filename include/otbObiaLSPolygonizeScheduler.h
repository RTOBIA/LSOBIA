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
#ifndef otbObiaLSPolygonizeScheduler_h
#define otbObiaLSPolygonizeScheduler_h

#include "otbObiaLSGraphToVectorScheduler.h"
#include "otbObiaGraphToVectorFilter.h"
#include "otbObiaSimplifyVectorFilter.h"
#include "otbOGRLayerWrapper.h"

/**
\file otbObiaLSPolygonizeScheduler.h
\brief This file define the class which will be used for filter processing a graph into a vector
 	 	in a multiprocess environnment.
 	 	Templated by:
 	 	- An input graph type
 	 	-A simplification function
*/
namespace otb
{
namespace obia
{

/**\class LSPolygonizeScheduler otbObiaLSPolygonizeScheduler.h
 * \brief Class handling multi-processing of a graph into a vector (simplified or not)*/
template< class TInputGraph, class TSimplifyFunc >
  class LSPolygonizeScheduler : public LSGraphToVectorScheduler< TInputGraph >
{
public:

    /** Standard class alias */
    using Self          = LSPolygonizeScheduler;
    using SuperClass   = itk::Object;
    using Pointer      = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >;

    /** Some convenient class alias */
	using InputGraphType   = TInputGraph;
  	using SimplifyFuncType = TSimplifyFunc;
	using SimplifyFunctionType = TSimplifyFunc;
	using InputGraphPointerType = typename InputGraphType::Pointer;
	using GraphOperationsType = GraphOperations<InputGraphType>;
	using OGRDataSourceType = otb::ogr::DataSource;
	using GraphToVectorFilterType = GraphToVectorFilter<InputGraphType>;
	using SimplifyFilterType = SimplifyVectorFilter<SimplifyFunctionType>;
	using OGRLayerType	 = otb::ogr::Layer;
	using OGRFeatureType = otb::ogr::Feature;
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(LSPolygonizeScheduler, itk::LightObject);

    /** Get/Set methods */
    void SetAggregateGraphs(bool AggregateGraphs){m_AggregateGraphs = AggregateGraphs;};

    //Set simplify
    void SetSimplifyFunc(SimplifyFuncType* simplifyFunc){ m_SimplifyFunc = simplifyFunc;};

    //Set simplify
    void SetIsSimplify(bool isSimplify){ m_IsSimplify = isSimplify;};
protected:

    /** Constructor */
    LSPolygonizeScheduler();

    /** Destructor */
    ~LSPolygonizeScheduler();

    /** Generation method */
    virtual void GenerateData();

    virtual void NoTilingExecution();

    virtual void TilingExecution();

private:

    /**\brief Create all filters required to convert graph into a fvector and simplify it if activated */
    void RunFilters();

    /**\brief Extract stability margins*/
    void ExtractStabilityMargins();

    /**\brief Aggregate stability margins to the graph*/
    void AggregateStabilityMargins();

    /**\brief Add metadata to the OGRDataset like tile width, tile height, image origin, etc...
     * \param: Tile giving metadata*/
    void AddMetaData(const ProcessingTile tile);

    /**\brief Remove all polygons outside the given tile. It checks among all polygons which one do not intersect
     * tile polygon and remove it from the layer m_OutputlayerName
     * \param: Tile to check*/
    void RemovePolygonsOutsideTile(const ProcessingTile& tile);

    /**\brief Create polygon version of a tile. It create a simple rectangular polygon using boundaries of the tile
     * \param: Tile
     * \return: A polygon representing the tile*/
    OGRPolygon* CreateTilePolygon(const ProcessingTile& tile);

    /**\brief Write features associated to the given tile. Name is set to Reconstructed_Polygons_ty_tx.gml
     * \param: Tile*/
    void WriteFeatures(const ProcessingTile& tile);

    /**\brief Convert graph to image, used for DEBUG
     * \param: Input graph
     * \param: Output filename*/
    void ConvertGraphToImage(InputGraphPointerType inputGraph, std::string filename);

    // The serialized stability margin extracted from the current graph
    std::vector< char > m_SerializedStabilityMargin;

    // The maximum number of bytes to send via MPI
    unsigned long int m_MaxNumberOfBytes;

    // Flag to activate or deactivate the aggreagation of graph
    bool m_AggregateGraphs;

    //Flag to indicate if we simplify vector
    bool m_IsSimplify;

    /**Output layer name*/
    std::string m_OutputLayerName;

    /**Current tile*/
    ProcessingTile m_CurrentTile;

    //Simplify func
    SimplifyFuncType* m_SimplifyFunc;

    //Tile processing
    bool m_IsTileProcessing;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaLSPolygonizeScheduler.txx"
#endif
