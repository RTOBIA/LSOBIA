#ifndef otbObiaLSGraphToGraphFilter_h
#define otbObiaLSGraphToGraphFilter_h

#include <cstdint>
#include <sstream>
#include "itkLightObject.h"
#include "otbObiaGraphOperations.h"


namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputGraph>
class LSGraphToGraphFilter : public itk::LightObject
{

public:

   /** Standard class alias */
   using Self 	      = LSGraphToGraphFilter;
   using SuperClass   = itk::Object;
   using Pointer      = itk::SmartPointer< Self >;
   using ConstPointer = itk::SmartPointer< const Self >;

   /** Run-time type information (and related methods). */
   itkTypeMacro(LSGraphToGraphFilter, itk::LightObject);

   /** Some convenient class alias */
   using InputGraphType = TInputGraph;
   using InputGraphPointerType = typename InputGraphType::Pointer;

   using OutputGraphType = TOutputGraph;
   using OutputGraphPointerType = typename OutputGraphType::Pointer;
   using GraphOperationsType = GraphOperations<OutputGraphType>;

   virtual void Update();


   /** Get/Set methods */
   //itkSetMacro(OutputDir, std::string);
   void SetAvailableMemory(const uint64_t mem);
   void SetOutputDir(const std::string path);
   void SetTemporaryDirectory(const std::string& tmpDir);
   void SetNumberOfTilesX(const uint32_t NumberOfTilesX){m_NumberOfTilesX = NumberOfTilesX;};
   void SetNumberOfTilesY(const uint32_t NumberOfTilesY){m_NumberOfTilesY = NumberOfTilesY;};
   void SetMaxTileSizeX(const uint32_t MaxTileSizeX){m_MaxTileSizeX = MaxTileSizeX;};
   void SetMaxTileSizeY(const uint32_t MaxTileSizeY){m_MaxTileSizeY = MaxTileSizeY;};
   void SetProjectionRef(const std::string ProjectionRef){m_ProjectionRef = ProjectionRef;};
   void SetImageWidth(const uint32_t ImageWidth){m_ImageWidth = ImageWidth;};
   void SetImageHeight(const uint32_t ImageHeight){m_ImageHeight = ImageHeight;};
   void SetTileMap(const std::map< uint32_t, ProcessingTile > tileMap){ m_TileMap = tileMap;};
   void SetTilesPerProcessor(const std::map< int, std::set<uint32_t> > tilesPerProcessor){ m_TilesPerProcessor = tilesPerProcessor;};
   void SetNumberOfSpectralBands(const uint32_t NumberOfSpectralBands){m_NumberOfSpectralBands = NumberOfSpectralBands;};
   void SetMaxNumberOfTilesPerProcessor(const uint32_t MaxNumberOfTilesPerProcessor){m_MaxNumberOfTilesPerProcessor = MaxNumberOfTilesPerProcessor;};

   /**TODO : Il faudra surement dupliquer le graph,
    *  car lorsque le filtre ayant généré ce graphe sera détruit, son graphe le sera aussi...*/
   void SetInputGraph(const InputGraphPointerType inputGraph){ m_InputGraph = inputGraph;};

 protected:

   LSGraphToGraphFilter();
   virtual ~LSGraphToGraphFilter();

   /**
	Abstract method that computes the initial width of padding crown to consider around
	the tiles to be processed.
    */
   virtual void ComputePaddingValue() = 0;

   /** Main method that needs to be override by the herited class */
   virtual void GenerateData() = 0;

   /** Output method that needs to be override by the herited class */
   virtual void CreateOutput() = 0;

   /** Compute all required parameters for processing*/
   virtual void PreProcessing();

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


   /** Display function*/
   void DisplayAttributes();

	/** Output directory where the final(s) graphs and the output label image (if necessary) are stored. */
	std::string m_OutputDir;
	bool m_SplittedGraph;
	bool m_WriteGraph;


	/** Available memory in bytes (in the master node) */
	uint64_t m_AvailableMemory;

	/** The temporary directory to store intermediate results */
	std::string m_TemporaryDirectory;

    /** The projection ref */
    std::string m_ProjectionRef;

    /** Number of tiles wrt to x-axis */
    uint32_t m_NumberOfTilesX;

    /** Number of tiles wrt to Y axis */
    uint32_t m_NumberOfTilesY;

    /** Size of tile X*/
    uint32_t m_MaxTileSizeX;

    /** Size of tile Y*/
    uint32_t m_MaxTileSizeY;

    /** Width of the input large scale image */
    uint32_t m_ImageWidth;

    /** Height of the input large scale image */
    uint32_t m_ImageHeight;

    /** Number of spectral bands of the input large scale image */
    uint32_t m_NumberOfSpectralBands;

	/** Tiles processed by this processor */
	std::map< uint32_t, ProcessingTile > m_TileMap;

	/** List of tiles per processor */
	std::map< int, std::set<uint32_t> > m_TilesPerProcessor;

	/** Maximum number of tiles for one processor */
	uint32_t m_MaxNumberOfTilesPerProcessor;

	/** Width of the crown of pixels to consider around the tiles */
	uint32_t m_PaddingValue;

	/**Input graph*/
	InputGraphPointerType m_InputGraph;

	/**Output Graph*/
	OutputGraphPointerType m_Graph;

private:


};
}//End obia
}//End otb

#include "otbObiaLSGraphToGraphFilter.txx"
#endif
