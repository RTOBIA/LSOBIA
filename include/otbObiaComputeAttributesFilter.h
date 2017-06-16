#ifndef otbObiaComputeAttributesFilter_h
#define otbObiaComputeAttributesFilter_h

#include <otbObiaGenericAttribute.h>
#include <otbPersistentSamplingFilterBase.h>

/**
\file otbObiaComputeAttributsFilter.h
\brief This file define the filter used to compute attributs. It extracts samples from a ROI and compute
all the attributes given
*/
namespace otb
{
namespace obia
{

/**Class specializing the merging cost function required by the generic filter*/
template <class TInputImage, class TMaskImage = otb::Image<unsigned char,2> >
class ComputeAttributesFilter : public PersistentSamplingFilterBase<TInputImage, TMaskImage>
{
	public:

		/** Standard typedefs */
	  using Self 				= ComputeAttributesFilter<TInputImage, TMaskImage>;
	  using Superclass   		= PersistentSamplingFilterBase<TInputImage, TMaskImage>;
	  using Pointer 			= itk::SmartPointer<Self>;
	  using ConstPointer	    = itk::SmartPointer<const Self>;
	  using OGRDataSourceType	= otb::ogr::DataSource;
	  using OGRDataPointer      = OGRDataSourceType::Pointer;
	  using OGRLayerType		= ogr::Layer;
	  using OGRFeatureType		= ogr::Feature;

	  /** Creation through object factory macro */
	  itkNewMacro(Self);

	  itkTypeMacro(ComputeAttributesFilter, PersistentSamplingFilterBase);


	  /** Template parameters typedefs */
	  using InputImageType 	     = TInputImage;
	  using MaskImageType 	     = TMaskImage;
	  using InputPixelType	     = typename InputImageType::InternalPixelType;
	  using RegionType 	  	     = typename InputImageType::RegionType;
	  using GenericAttributeType = GenericAttribute<InputImageType>;


	  /** Set/Get macro for the field name containing class names
	    * in the input vectors.*/
	   itkSetMacro(InputLayerName, std::string);
	   itkGetMacro(InputLayerName, std::string);

	   itkSetMacro(OutputDir, std::string);
	   itkGetMacro(OutputDir, std::string);

	   itkSetMacro(OutputFilename, std::string);
	   itkGetMacro(OutputFilename, std::string);

	   /**
	    * Reset the persistent data of the filter.
	    */
	   virtual void Reset(void);
	   /**
	    * Synthesize the persistent data of the filter.
	    */
	   virtual void Synthetize(void);

	  /**Set all computing attributs method*/
	   void SetAttributs(std::vector<GenericAttributeType*> attributs){ m_Attributs = attributs;};

	protected:
	  /** Constructor */
	  ComputeAttributesFilter();
	  /** Destructor */
	  ~ComputeAttributesFilter() ITK_OVERRIDE {};

	  /** Generate data should thread over */
	  void GenerateData(void) ITK_OVERRIDE;

	  /** Compute all attributs*/
	  void ComputeAllAttributs();

	  /**Initialize output */
	  void InitializeOutput();

	  /** Generic method called for each matching pixel position*/
	  virtual void ProcessSample(const ogr::Feature& feature,
	                             typename TInputImage::IndexType& imgIndex,
	                             typename TInputImage::PointType& imgPoint,
	                             itk::ThreadIdType& threadid) ITK_OVERRIDE;

	  /**Update the layer **/
	  void UpdateLayer(OGRLayer* layer, ogr::Feature& feature,
			  	  	   std::string fieldName, double fieldValue);

	private:

	  //Input layer name used for computing attributs
	  std::string m_InputLayerName;

	  //Vector of GenericAttributs
	  std::vector<GenericAttributeType*> m_Attributs;

	  //Vector containing samples
	  std::map<unsigned int, std::vector<InputPixelType> > m_Samples;

	  //Input image
	  InputImageType* m_InputImage;

	  //Number of components
	  unsigned int m_NumberOfBands;

	  //Output DS
	  OGRDataSourceType::Pointer m_OutputDs;

	  //Output dir
	  std::string m_OutputDir;

	  //Output file name
	  std::string m_OutputFilename;



};

} // end of namespace obia
} // end of namespace otb
#include <otbObiaComputeAttributesFilter.txx>
#endif

