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
