#ifndef otbObiaGraphFileReader_h
#define otbObiaGraphFileReader_h

#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkDataObject.h"

#include "otbObiaGraphSource.h"

#include <string>
#include <ostream>

namespace otb
{
namespace obia
{

/** \class GraphFileReaderException
 *
 * \brief Base exception class for IO conflicts.
 *
 */
class GraphFileReaderException : public itk::ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( GraphFileReaderException, ExceptionObject );

  /** Constructor. */
  GraphFileReaderException(const char *file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") :
    ExceptionObject(file, line, message, loc)
  {
  }

  /** Constructor. */
  GraphFileReaderException(const std::string &file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") :
    ExceptionObject(file, line, message, loc)
  {
  }
};


/** \class GraphFileReader
 *	\brief Base class for all process objects that output graph data.
 *
 * GraphFileReader is the base class for all process objects that output
 * graph data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the adjacent graph.
 */
template< typename TOutputGraph >
class GraphFileReader : public otb::obia::GraphSource<TOutputGraph>
{
public:

	/** Standard class alias. */
	using Self         = GraphFileReader;
	using Superclass   = itk::ProcessObject;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);
	/** Run-time type information (and related methods). */
	itkTypeMacro(GraphFileReader, ProcessObject);

	/** Some convenient alias */
	using OutputGraphType = TOutputGraph;

	/** */
	void SetFileName(std::string extendedFileName)
	{
		m_FileName = extendedFileName;
	}

protected:

    /** Does the real work */
	virtual void GenerateData();
	GraphFileReader();
  	virtual ~GraphFileReader() {}
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

	virtual void UpdateOutputInformation();
private:
    /** Test whether the given filename exist and it is readable,
        this is intended to be called before attempting to use
        ImageIO classes for actually reading the file. If the file
        doesn't exist or it is not readable, and exception with an
        appropriate message will be thrown. */
    void TestFileExistanceAndReadability();

    GraphFileReader(const Self &); //purposely not implemented
    void operator =(const Self&); //purposely not implemented

    std::string m_FileName; // The file to be read
};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphFileReader.txx"
#endif
