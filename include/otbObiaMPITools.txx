/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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
