#ifndef otbObiaGraphSource_h
#define otbObiaGraphSource_h
#include "itkProcessObject.h"
#include "otbObiaGraph.h"

namespace otb
{
namespace obia
{

/** \class GraphSource
 *	\brief Base class for all process objects that output graph data.
 *
 * GraphSource is the base class for all process objects that output
 * graph data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the adjacent graph.
 */
template< typename TOutputGraph >
class GraphSource : public itk::ProcessObject
{

public:

	/** Standard class alias. */
	using Self         = GraphSource;
	using Superclass   = itk::ProcessObject;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Smart Pointer type to a DataObject. */
	using DataObjectPointer = itk::DataObject::Pointer;

	using DataObjectIdentifierType       = Superclass::DataObjectIdentifierType;
	using DataObjectPointerArraySizeType = Superclass::DataObjectPointerArraySizeType;

	/** Run-time type information (and related methods). */
	itkTypeMacro(GraphSource, ProcessObject);

	/** Some convenient alias */
	using OutputGraphType = TOutputGraph;

	OutputGraphType * GetOutput();
	const OutputGraphType * GetOutput() const;
	OutputGraphType * GetOutput(unsigned int idx);
	OutputGraphType* GetOutputByMove();

	virtual itk::ProcessObject::DataObjectPointer 
			MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  virtual itk::ProcessObject::DataObjectPointer 
  		MakeOutput(const itk::ProcessObject::DataObjectIdentifierType &) ITK_OVERRIDE;

  virtual void GraftOutput(itk::DataObject *output);
  virtual void GraftOutput(const DataObjectIdentifierType & key, itk::DataObject *output);
  virtual void GraftNthOutput(unsigned int idx, itk::DataObject *output);

protected:

	GraphSource();
  	virtual ~GraphSource() {}

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphSource.txx"
#endif
