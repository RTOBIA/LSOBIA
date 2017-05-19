#ifndef otbObiaSimplifyVectorFilter_h
#define otbObiaSimplifyVectorFilter_h

#include "itkProcessObject.h"
#include <otbOGRDataSourceWrapper.h>

namespace otb
{
namespace obia
{

/** \class GraphToGraphFilter
 *	\brief Base class for all process objects that requires one input graph
 	and that provides one output graph data.
 *
 */
template< typename TSimplifyFunc>
class ITK_EXPORT SimplifyVectorFilter : public itk::ProcessObject
{

public:

	/** Standard class alias */
	using Self         = SimplifyVectorFilter;
	using Superclass   = itk::ProcessObject;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Run-time type information (and related methods). */
  	itkNewMacro(Self);

  	/** Some convenient typedefs. */
  	using SimplifyFunc = TSimplifyFunc;

  	/** Return the name of the class. */
	itkTypeMacro(SimplifyVectorFilter, ProcessObject);

	  /** Definition of the input image */
	using OGRDataSourceType = otb::ogr::DataSource;
	using OGRDataSourcePointerType = OGRDataSourceType::Pointer;
	using OGRLayerType = otb::ogr::Layer;
	using OGRFeatureType = otb::ogr::Feature;
	using DataObjectPointerArraySizeType =  itk::ProcessObject::DataObjectPointerArraySizeType;

	/** Set/Get the input graph of this process object.  */
	using Superclass::SetInput;
	virtual void SetInput(const OGRDataSourceType *input);
	virtual const OGRDataSourceType * GetInput(void);

	/** Set the Field Name in which labels will be written. (default is "DN")
	* A field "FieldName" of type integer is created in the output memory layer.
	*/
	itkSetMacro(FieldName, std::string);
	/**
	* Return the Field name in which labels have been written.
	*/
	itkGetMacro(FieldName, std::string);

	/**
	* Get the layer name in which features have been written.
	*/
	itkSetMacro(LayerName, std::string);

	/**
	* Get the output \c ogr::DataSource which is a "memory" datasource.
	*/
	const OGRDataSourceType * GetOutput();

	/**
	 * Set the simplify function*/
	void SetSimplifyFunc(SimplifyFunc* simplifyFunc){m_simplifyFunc = simplifyFunc;};

	protected:

	SimplifyVectorFilter();
	~SimplifyVectorFilter() ITK_OVERRIDE {}

	void GenerateInputRequestedRegion() ITK_OVERRIDE;

	/** Generate Data method*/
	void GenerateData() ITK_OVERRIDE;

	/** DataObject pointer */
	typedef itk::DataObject::Pointer DataObjectPointer;

	/**Create output*/
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

	using Superclass::MakeOutput;

	/**GDAL Method*/
	GDALDriver* initializeGDALDriver(std::string driverName);

	/**Initialize output DS*/
	void InitializeOutputDS();

	/**Create bounding feature to manage border*/
	OGRFeatureType* CreateBoundingBoxFeature();

	/**Intersect with bounding box*/
	std::vector<OGRGeometry*> IntersectWithBB(OGRFeatureType feature);

	/**Convert feature to edge
	 * @param : Current feature
	 * @param: adjacent features
	 * @param: Layer to store englobed polygons, will be processed later*/
	std::vector<OGRFeatureType*> ConvertToEdges(OGRFeatureType feature, std::vector<OGRFeatureType> adjFeatures, OGRLayerType& insidePolygons);

	/**Create an edge*/
	OGRFeatureType* CreateEdge(OGRGeometry* intersectedGeometry,
						   	   OGRFeatureType refFeature, OGRFeatureType adjFeature,
							   bool isSimplify = true);

	/**Compute interected features*/
	void IntersectFeatures(std::vector<OGRFeatureType*>& edges,
						   OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify = true);

	/**Compute Edge*/
	void ConvertGeometryCollectionToEdge(OGRGeometry* intersectedLine, std::vector<OGRFeatureType*>& edges,
										 OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify = true);

	/**Reconstruct all polygons from lines geometry*/
	void ReconstructAllPolygons();

	/** Reconstruct polygon*/
	OGRFeatureType* ReconstructPolygon(double startCoords);

	/** Sort linestring*/
	std::vector<OGRLineString*> ConvertToGeometries(double startCoords, std::vector<double>& adjCoords);

	//Write into a file
	void WriteFile(std::string filepath, int layerId = -1);

	private:

	SimplifyVectorFilter(const Self &);  //purposely not implemented

	void operator =(const Self&);      //purposely not implemented

	std::string m_FieldName;
	SimplifyFunc* m_simplifyFunc;
	std::string m_LayerName;

	OGRFeatureType* m_BBFeature;

	double m_CurrentCoords;
	//Edges map
	std::map<double, std::vector<OGRFeatureType*>> m_PolygonEdges;

	//BB Edges
	std::map<double, std::vector<OGRGeometry*>> m_bbEdges;


};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaSimplifyVectorFilter.txx"
#endif
