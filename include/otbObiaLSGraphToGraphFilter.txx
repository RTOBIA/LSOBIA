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
#ifndef otbObiaLSGraphToGraphFilter_txx
#define otbObiaLSGraphToGraphFilter_txx
#include "otbObiaLSGraphToGraphFilter.h"


namespace otb
{
namespace obia
{
template< typename TInputGraph, typename TOutputGraph>
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::LSGraphToGraphFilter() :
m_AvailableMemory(0),
m_TemporaryDirectory(""),
m_OutputDir(""),
m_SplittedGraph(false),
m_WriteGraph(false),
m_MaxNumberOfTilesPerProcessor(0),
m_PaddingValue(0),
m_InputGraph(nullptr),
m_Graph(nullptr)
{
}

template< typename TInputGraph, typename TOutputGraph>
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::~LSGraphToGraphFilter()
{
}

template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::PreProcessing()
{
    //Graft the graph
    this->m_Graph  = this->m_InputGraph;

    //Display attributes
    this->DisplayAttributes();

    //Compute Padding value
    this->ComputePaddingValue();
}

template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::Update()
{

    /**Pre-processing*/
    this->PreProcessing();

    /**Generate the data*/
    this->GenerateData();

    /**Create output*/
    this->CreateOutput();

}


template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::WriteGraphIfNecessary(const unsigned int ty,
                            const unsigned int tx)
{
    if(m_TileMap.size() > 1)
    {
        // There are other tiles to process, we need to store this one
        // on the local disk of this processor.
        std::stringstream os;
        os << m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
        GraphOperationsType::WriteGraphToDisk(m_Graph, os.str());
        m_Graph->Reset();
    }
}


template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::ReadGraphIfNecessary(const unsigned int ty,
                       const unsigned int tx)
{
    if(m_TileMap.size() > 1)
    {
        std::stringstream in;
        in << m_TemporaryDirectory << "Graph_" << ty << "_" << tx << ".dat";
        m_Graph = GraphOperationsType::ReadGraphFromDisk(in.str());
    }
}

template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::SetOutputDir(const std::string path)
{
    m_OutputDir = path;
    if( !m_OutputDir.empty() && m_OutputDir[m_OutputDir.size() - 1] != '/')
    {
        m_OutputDir.append("/");
    }
}


template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::SetAvailableMemory(const uint64_t mem)
{
    m_AvailableMemory = mem * 1024 * 1024;
}

template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::SetTemporaryDirectory(const std::string& tmpDir)
{
    m_TemporaryDirectory = tmpDir;

    if(!m_TemporaryDirectory.empty())
    {
        if(m_TemporaryDirectory[m_TemporaryDirectory.size()-1] != '/')
        {
            m_TemporaryDirectory.append("/");
        }
    }
}

template< typename TInputGraph, typename TOutputGraph>
void
LSGraphToGraphFilter<TInputGraph, TOutputGraph>
::DisplayAttributes()
{
    std::cout << " ------- Class Attributes : --------" << std::endl;
    std::cout << "Available Memory : " << m_AvailableMemory << std::endl;
    std::cout << "Temporary Directory : " << m_TemporaryDirectory << std::endl;
    std::cout << "Projection Ref : " << m_ProjectionRef << std::endl;
    std::cout << "Number of Tiles X : " << m_NumberOfTilesX << std::endl;
    std::cout << "Number of Tiles Y: " << m_NumberOfTilesY << std::endl;
    std::cout << "Tile Size X : " << m_MaxTileSizeX << std::endl;
    std::cout << "Tile Size Y : " << m_MaxTileSizeY << std::endl;
    std::cout << "Image Width : " << m_ImageWidth << std::endl;
    std::cout << "Image Height : " << m_ImageHeight << std::endl;
    std::cout << "Number of Spectral Bands : " << m_NumberOfSpectralBands << std::endl;
    std::cout << "Max Number Of Tiles : " << m_MaxNumberOfTilesPerProcessor << std::endl;
    if(m_Graph != nullptr)
    {
        std::cout << "Nombre de noeuds du graphe : " << m_Graph->GetNumberOfNodes() << std::endl;
    }

}

}//end obia

}//end otb
#endif
