#ifndef otbObiaLSPolygonizeScheduler_h
#define otbObiaLSPolygonizeScheduler_h

#include "otbObiaLSGraphToVectorScheduler.h"
#include "otbObiaGraphToVectorFilter.h"
#include "otbObiaSimplifyVectorFilter.h"
#include "otbOGRLayerWrapper.h"

namespace otb
{
namespace obia
{

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

    void RunFilters();
    void RescaleGraph(ProcessingTile& tile);
    void ExtractStabilityMargins();
    void AggregateStabilityMargins();
    void RemovePolygonsOutsideTile(ProcessingTile& tile);
    OGRPolygon* CreateTilePolygon(ProcessingTile& tile);

    //Write the features belonging to a tile
    void WriteFeatures(ProcessingTile& tile);

    //For debug
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

    //Simplify func
    SimplifyFuncType* m_SimplifyFunc;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaLSPolygonizeScheduler.txx"
#endif
