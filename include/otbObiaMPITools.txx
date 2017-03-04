#ifndef otbObiaMPITools_txx
#define otbObiaMPITools_txx
#include "otbObiaMPITools.h"

namespace otb
{
namespace obia
{
/** Initialize the singleton */
MPITools::Pointer MPITools::m_Singleton = nullptr;


MPITools::Pointer MPITools::Instance()
{
  if(m_Singleton.GetPointer() == NULL)
    {
    m_Singleton = itk::ObjectFactory<Self>::Create();

    if(m_Singleton.GetPointer() == NULL)
      {
      m_Singleton = new MPITools;
      }
    m_Singleton->UnRegister();
    }

  return m_Singleton;
}

MPITools::MPITools()
{
}

MPITools::~MPITools()
{
}

bool 
MPITools::IsMyTurn(const uint32_t id)
{
	auto mpiConfig = otb::MPIConfig::Instance();
	return ( mpiConfig->GetMyRank() == (id % mpiConfig->GetNbProcs()) );
}

template< typename T >
void 
MPITools::Accumulate(T& accvalue, MPI_Datatype datatype)
{
	auto mpiConfig = otb::MPIConfig::Instance();

	// Master process
	if(mpiConfig->GetMyRank() == 0)
	{
		// Receive from other processes the number of bytes
		for(uint32_t p = 1; p < mpiConfig->GetNbProcs(); p++)
		{
			T value;
			MPI_Recv(&value, 1, datatype, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			accvalue += value;
		}

		// Send accumulated value to all slaves
		for(uint32_t p = 1; p < mpiConfig->GetNbProcs(); p++)
		{
			MPI_Send(&accvalue, 1, datatype, p, 0, MPI_COMM_WORLD);
		}

	}
	// Slave processes
	else
	{
		MPI_Send(&accvalue, 1, datatype, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(&accvalue, 1, datatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

int 
MPITools::GetProcessorRankFromTileId(const uint32_t tid)
{
	auto mpiConfig = otb::MPIConfig::Instance();
	return (tid % mpiConfig->GetNbProcs());
}

template< typename T >
void 
MPITools::ComputeMax(T& maxValue, MPI_Datatype datatype)
{
	auto mpiConfig = otb::MPIConfig::Instance();

	// Master process
	if(mpiConfig->GetMyRank() == 0)
	{
		// Receive from other processes the number of bytes
		for(uint32_t p = 1; p < mpiConfig->GetNbProcs(); p++)
		{
			T value;
			MPI_Recv(&value, 1, datatype, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			maxValue = std::max(value, maxValue);
		}

		// Dispatch
		for(uint32_t p = 1; p < mpiConfig->GetNbProcs(); p++)
		{
			MPI_Send(&maxValue, 1, datatype, p, 0, MPI_COMM_WORLD);
		}
	}
	// Slave processes
	else
	{
		MPI_Send(&maxValue, 1, datatype, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(&maxValue, 1, datatype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}


} // end of namespace obia
} // end of namespace otb

#endif
