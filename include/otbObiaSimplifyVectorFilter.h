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

	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;
	using Superclass::MakeOutput;

	/**GDAL Method*/
	GDALDriver* initializeGDALDriver(std::string driverName);


	//Write into a file
	void WriteFile(std::string filepath);

	private:

	SimplifyVectorFilter(const Self &);  //purposely not implemented
	void operator =(const Self&);      //purposely not implemented

	std::string m_FieldName;
	SimplifyFunc* m_simplifyFunc;

};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaSimplifyVectorFilter.txx"
#endif
