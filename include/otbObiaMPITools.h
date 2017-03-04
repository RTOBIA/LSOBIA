#ifndef otbObiaMPITools_h
#define otbObiaMPITools_h
#include <cstdint>
#include "itkObject.h"
#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "mpi.h"
#include "otbMPIConfig.h"

namespace otb
{
namespace obia
{

class MPITools : public itk::LightObject
{
public:

	/** Standard class alias */
	using Self         = MPITools;
	using Superclass   = itk::Object;
	using Pointer      = itk::SmartPointer< Self >;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Retrieve the singleton instance */
	static Pointer Instance();

	/** Run-time type information (and related methods). */
	itkTypeMacro(MPITools, itk::LightObject);

	/** Given the rank, the tile (or graph) index and the number of processors
        this method indicate if the current processor is in charge of the current tile */
	bool IsMyTurn(const uint32_t id);

	/** Compute the sum of all the values across the cluster and send the
 		sum to the slaves */
	template< typename T >
	void Accumulate(T& accvalue, MPI_Datatype datatype);

	/** Compute the max of all the values across the cluster and send it
		to the slaves */
	template< typename T >
	void ComputeMax(T& maxValue, MPI_Datatype datatype);

	/** Given the tile id, this methods returns the matched processor. */
	int GetProcessorRankFromTileId(const uint32_t tid);

protected:

	/** Constructor */
	MPITools();

	/** Destructor */
	virtual ~MPITools();

private:

  MPITools(const MPITools &); //purposely not implemented
  void operator =(const MPITools&); //purposely not implemented


  static Pointer m_Singleton;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaMPITools.txx"
#endif
