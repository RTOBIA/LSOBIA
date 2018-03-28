#ifndef otbObiaLSImageToGraphScheduler_h
#define otbObiaLSImageToGraphScheduler_h
#include <cstdint>
#include <sstream>

#include "otbObiaGraphOperations.h"
#include "otbObiaMPITools.h"
#include "otbObiaGraphToLabelImageFilter.h"
#include "otbImageFileReader.h"
#include "otbMultiChannelExtractROI.h"
#include "otbMPIConfig.h"
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImage.h"
#include "otbImageFileWriter.h"


/**
\file otbObiaLSImageToGraphScheduler.h
\brief This file define the class which will be used for filter processing an image to a graph
 	 	in a multiprocess environnment
*/
namespace otb
{
namespace obia
{

/**\class LSImageToGraphScheduler otbObiaLSImageToGraphScheduler.h
 * \brief Class handling multi-processing of an image to a graph
 * Templated by:
 * - an input image type
 * - an output graph type (like baatz)*/
template< typename TInputImage, typename TOutputGraph>
class LSImageToGraphScheduler : public itk::LightObject
{
public:
    
    /** Standard class alias */
    using Self         = LSImageToGraphScheduler;
    using SuperClass   = itk::Object;
    using Pointer      = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >; 

    /** Run-time type information (and related methods). */
    itkTypeMacro(LSImageToGraphScheduler, itk::LightObject);

    /** Some convenient class alias */
    using InputImageType = TInputImage;
    using InputImageReaderType = ImageFileReader<InputImageType>;
    using OutputGraphType = TOutputGraph;
    using OutputGraphPointerType = typename OutputGraphType::Pointer;
    using GraphOperationsType = GraphOperations<OutputGraphType>;

    virtual void Update();


    /** Get/Set methods */
    void SetFileName(const std::string& filename);
    void SetWriteLabelImage(bool v){ m_WriteLabelImage = v;}
    void SetWriteGraph(bool v){ m_WriteGraph = v;}
    void SetOutputDir(const std::string & path);
    void SetMaxTileSizeX(const uint32_t MaxTileSizeX);
    void SetMaxTileSizeY(const uint32_t MaxTileSizeY);
    void SetLabelImageName(const std::string& labelImage);
    void SetProcessNoData(bool p){ m_ProcessNoData = p;}
    void SetNoDataValue(const float n);


    // The memory is given in megabytes.
    void SetAvailableMemory(const uint64_t mem);
    void SetTemporaryDirectory(const std::string& tmpDir);

    /**TODO Get the output graph : maybe duplicate the graph*/
    OutputGraphPointerType GetGraph(){return m_Graph;};
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
    bool GetProcessNoData(){return m_ProcessNoData;};
    bool GetNoDataValue(){return m_NoDataValue;};

protected:

    LSImageToGraphScheduler();
    virtual ~LSImageToGraphScheduler();

    /** 
    Abstract method that computes the initial width of padding crown to consider around
    the tiles to be processed.
     */
    virtual void ComputePaddingValue() = 0;

    /** Main method that needs to be override by the herited class */
    virtual void GenerateData() = 0;

    /** Prepare the tiles to be processed. */
    void PreProcessing();

    /** 
       Simple method that writes a graph if the number of tiles to be processed per
       processor is greater than 1.
     */
    void WriteGraphIfNecessary(const unsigned int ty, 
                               const unsigned int tx);

    /**
        Simple method that loads a graph if the number of tiles per processor is
        greater than 1.
    */
    void ReadGraphIfNecessary(const unsigned int ty,
                              const unsigned int tx);

    /**Write graph into a label image*/
    void ConvertGraphToImage(const unsigned int ty,
            				 const unsigned int tx);


    /** Path of the input big image */
    std::string m_FileName;

    /** Output directory where the final(s) graphs and the output label image (if necessary) are stored. */
    std::string m_OutputDir;
    std::string m_LabelImageName;
    bool m_WriteLabelImage;
    bool m_SplittedGraph;
    bool m_WriteGraph;
    bool m_ProcessNoData;

    /** Tile size with respect to x-axis. */
    uint32_t m_MaxTileSizeX;
    
    /** Tile size with respect to y-axis. */
    uint32_t m_MaxTileSizeY;

    /** Available memory in bytes (in the master node) */
    uint64_t m_AvailableMemory;

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

    /**Image origin X*/
    uint32_t m_OriginX;

    /**Image origin Y*/
    uint32_t m_OriginY;

    /** Tiles processed by this processor */
    std::map< uint32_t, ProcessingTile > m_TileMap;
    
    /** List of tiles per processor */
    std::map< int, std::set<uint32_t> > m_TilesPerProcessor;
    
    /** Maximum number of tiles for one processor */
    uint32_t m_MaxNumberOfTilesPerProcessor;

    /** Width of the crown of pixels to consider around the tiles */
    uint32_t m_PaddingValue;

    OutputGraphPointerType m_Graph;

    /** Value set for NO DATA */
    float m_NoDataValue;

private:
    
    /** Compute the final tile size by optimizing load balancing */
    void ComputeFinalTileSize();

    /** Assign tiles to each processor and compute the tile frames */
    void ComputeTileFrames();

    /** Create all necessary files to the output directory */
    void CreateOutput();
    
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaLSImageToGraphScheduler.txx"
#endif
