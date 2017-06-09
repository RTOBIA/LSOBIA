#ifndef otbObiaLSGraphToVectorScheduler_h
#define otbObiaLSGraphToVectorScheduler_h
#include <cstdint>
#include <sstream>

#include "otbObiaGraphOperations.h"
#include <otbOGRDataSourceWrapper.h>

namespace otb
{
namespace obia
{

template< typename TInputGraph>
class LSGraphToVectorScheduler : public itk::LightObject
{
public:

    /** Standard class alias */
    using Self         = LSGraphToVectorScheduler;
    using SuperClass   = itk::Object;
    using Pointer      = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >;

    /** Run-time type information (and related methods). */
    itkTypeMacro(LSGraphToVectorScheduler, itk::LightObject);

    /** Some convenient class alias */
    using InputGraphType = TInputGraph;
    using InputGraphPointerType = typename InputGraphType::Pointer;
    using OGRLayerType			= typename otb::ogr::Layer;
    using GraphOperationsType = GraphOperations<InputGraphType>;
    using OGRDataSourceType = otb::ogr::DataSource;
    using OGRDataSourcePointerType = otb::ogr::DataSource::Pointer;

    virtual void Update();


    /** Get/Set methods */
    void SetFileName(const std::string& filename);
    void SetWriteVector(bool v){ m_WriteVector = v;}
    void SetWriteGraph(bool v){ m_WriteGraph = v;}
    void SetOutputDir(const std::string & path);
    void SetMaxTileSizeX(const uint32_t MaxTileSizeX);
    void SetMaxTileSizeY(const uint32_t MaxTileSizeY);
    void SetGraphPrefixName(const std::string prefixName){m_GraphPrefixName = prefixName;};

    // The memory is given in megabytes.
    void SetAvailableMemory(const uint64_t mem);
    void SetInputDirectory(const std::string& inputDir);
    void SetTemporaryDirectory(const std::string& tmpDir);

    /**TODO Get the output graph : maybe duplicate the graph*/
    InputGraphPointerType GetGraph(){return m_Graph;};
    std::map< uint32_t, ProcessingTile > GetTileMap(){return m_TileMap;};
    std::map< int, std::set<uint32_t> > GetTilesPerProcessor(){return m_TilesPerProcessor;};
    uint32_t GetImageWidth(){return m_ImageWidth;};
    uint32_t GetImageHeight(){return m_ImageHeight;};
    std::string GetProjectionRef(){return m_ProjectionRef;};
    uint32_t GetNumberOfTilesX(){return m_NumberOfTilesX;};
    uint32_t GetNumberOfTilesY(){return m_NumberOfTilesY;};
    uint32_t GetMaxTileSizeX(){return m_MaxTileSizeX;};
    uint32_t GetMaxTileSizeY(){return m_MaxTileSizeY;};
    uint32_t GetNumberOfSpectralBands(){return m_NumberOfSpectralBands;};
    uint32_t GetMaxNumberOfTilesPerProcessor(){return m_MaxNumberOfTilesPerProcessor;};
    void SetTileMap(const std::map< uint32_t, ProcessingTile > & tileMap){ m_TileMap = tileMap;};
    void SetTilesPerProcessor(std::map< int, std::set<uint32_t> > tilesPerProcessor){ m_TilesPerProcessor = tilesPerProcessor;};
    void SetGraph(InputGraphPointerType graph){m_Graph = graph;};
protected:

    LSGraphToVectorScheduler();
    virtual ~LSGraphToVectorScheduler();

    /** Main method that needs to be override by the herited class */
    virtual void GenerateData() = 0;

    /** Prepare the tiles to be processed. */
    void PreProcessing();

    /** Compute tiles frame*/
    void ComputeTileFrames();
    /**
       Simple method that writes a graph if the number of tiles to be processed per
       processor is greater than 1.
     */
    void WriteGraphIfNecessary(const unsigned int ty,
                               const unsigned int tx,
							   std::string graph_dir,
							   std::string prefix_name);

    /**
        Simple method that loads a graph if the number of tiles per processor is
        greater than 1.
    */
    void ReadGraphIfNecessary(const unsigned int ty,
                              const unsigned int tx,
							  std::string graph_dir,
							  std::string prefix_name);

    /**
     * Simple method that writes a vector into a gml files that can be used after*/
    void WriteVectorIfNecessary(const unsigned int ty,
                               const unsigned int tx);

    /**
         Simple method that read a vector if the number of tiles per processor is
         greater than 1.
     */
     void ReadVectorIfNecessary(const unsigned int ty,
                               const unsigned int tx);

    /** Path of the input big image */
    std::string m_FileName;

    /** Output directory where the final(s) graphs and the output vector (if necessary) are stored. */
    std::string m_OutputDir;
    std::string m_VectorName;
    bool m_WriteVector;
    bool m_SplittedGraph;
    bool m_WriteGraph;

    /** Tile size with respect to x-axis. */
    uint32_t m_MaxTileSizeX;

    /** Tile size with respect to y-axis. */
    uint32_t m_MaxTileSizeY;

    /** Available memory in bytes (in the master node) */
    uint64_t m_AvailableMemory;

    /** Input directory where graph are stored*/
    std::string m_InputDirectory;

    /** The temporary directory to store intermediate results */
    std::string m_TemporaryDirectory;

    /** Width of the input large scale image */
    uint32_t m_ImageWidth;

    /** Height of the input large scale image */
    uint32_t m_ImageHeight;

    /** Number of spectral bands of the input large scale image */
    uint32_t m_NumberOfSpectralBands;

    /** The projection ref */
    std::string m_ProjectionRef;

    /** Number of tiles wrt to x-axis */
    uint32_t m_NumberOfTilesX;

    /** Number of tiles wrt to Y axis */
    uint32_t m_NumberOfTilesY;

    /** Tiles processed by this processor */
    std::map< uint32_t, ProcessingTile > m_TileMap;

    /** List of tiles per processor */
    std::map< int, std::set<uint32_t> > m_TilesPerProcessor;

    /** Maximum number of tiles for one processor */
    uint32_t m_MaxNumberOfTilesPerProcessor;

    //Padding value
    unsigned int m_PaddingValue;

    /** Graph prefix name*/
    std::string m_GraphPrefixName;

    //Input graph
    InputGraphPointerType m_Graph;

    //Output OGR DS
    OGRDataSourcePointerType m_OutputDS;

private:

    /** Create all necessary files to the output directory */
    void CreateOutput();

};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaLSGraphToVectorScheduler.txx"
#endif
