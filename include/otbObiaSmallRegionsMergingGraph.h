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
#ifndef otbObiaSmallRegionsMergingGraph_h
#define otbObiaSmallRegionsMergingGraph_h
#include <limits>

#include "otbObiaGraph.h"

namespace otb
{
namespace obia
{

struct SRMEdgeAttribute : GraphAttribute
{
    // Merging cost
    float m_MergingCost;

    // Flag indicating if the merging cost has to be recomputed.
    bool m_CostUpdated;

    SRMEdgeAttribute(){}

    SRMEdgeAttribute(const SRMEdgeAttribute& other)
        : m_MergingCost(other.m_MergingCost), m_CostUpdated(other.m_CostUpdated)
    {}

    virtual uint64_t GetMemorySize() const;

    virtual uint64_t GetNumberOfBytesToSerialize() const;

    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const
    {
        AvoidCompilationWarings(stream, position);
    }

    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position)
    {
        AvoidCompilationWarings(stream, position);
    }

private:
    // Does nothing
    virtual void AvoidCompilationWarings(const std::vector<char>& stream, uint64_t& position) const
    {
        // Avoid compilation warnings...
        (void)stream;
        (void)position;
    }

};

struct SRMNodeAttribute : GraphAttribute
{
    // Average spectral values for each band
    std::vector< float > m_AvgSpec;

    // Standard deviation for each spectral band
    std::vector< float > m_StdSpec;

    // Area of the segment
    uint32_t m_Area;

    // Perimeter of the segment
    uint32_t m_Perimeter;

    // To indicate if the node has merged at the previous iteration
    bool m_HasPreviouslyMerged;

    SRMNodeAttribute(){}

    SRMNodeAttribute(const SRMNodeAttribute& other)
    : m_AvgSpec(other.m_AvgSpec), m_StdSpec(other.m_StdSpec),
      m_Area(other.m_Area), m_Perimeter(other.m_Perimeter),
      m_HasPreviouslyMerged(other.m_HasPreviouslyMerged)
    {}

    virtual uint64_t GetMemorySize() const;

    virtual uint64_t GetNumberOfBytesToSerialize() const;

    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const;

    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position);
};


}//End otb
}//End obia
#include "otbObiaSmallRegionsMergingGraph.txx"
#endif
