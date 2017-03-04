
#ifndef otbObiaSmallRegionsMergingGraph_txx
#define otbObiaSmallRegionsMergingGraph_txx
#include "otbObiaSmallRegionsMergingGraph.h"
namespace otb
{
namespace obia
{


uint64_t
SRMEdgeAttribute::GetMemorySize() const
{
	return FloatSize + BoolSize;
}

uint64_t
SRMEdgeAttribute::GetNumberOfBytesToSerialize() const
{
	// No specific data needs to be serialized
	return 0;
}

uint64_t
SRMNodeAttribute::GetMemorySize() const
{
  return m_AvgSpec.size() * 2 * FloatSize + // Mean and standard deviations
  	  2 * sizeof(std::vector<float>) +
	  2 * UInt32Size + // area and perimeter
	  BoolSize; // the flag
}

uint64_t
SRMNodeAttribute::GetNumberOfBytesToSerialize() const
{
  // The flag does not need to be serialized.
  return UInt32Size +
    m_AvgSpec.size() * 2 * FloatSize +
    2* UInt32Size;
}

void
SRMNodeAttribute::Serialize(std::vector<char>& stream, uint64_t& position) const
{
	// Serialize the number of bands
	const uint32_t numBands = m_AvgSpec.size();
	std::memcpy(&(stream[position]), &numBands, UInt32Size);
	position += UInt32Size;

	// Serialize the average spectral values
	std::memcpy(&(stream[position]), &m_AvgSpec[0], numBands * FloatSize);
	position += numBands * FloatSize;

	// Serialize the standard deviation values
	std::memcpy(&(stream[position]), &m_StdSpec[0], numBands * FloatSize);
	position += numBands * FloatSize;

	// Serialize the perimeter
	std::memcpy(&(stream[position]), &m_Perimeter, UInt32Size);
	position += UInt32Size;

	// Serialize the area
	std::memcpy(&(stream[position]), &m_Area, UInt32Size);
	position += UInt32Size;
}

void
SRMNodeAttribute::DeSerialize(const std::vector<char>& stream, uint64_t& position)
{
	// Deserialize the number of bands.
	uint32_t numBands = 0;
	std::memcpy(&numBands, &(stream[position]), UInt32Size);
	position += UInt32Size;

	// Initialization of the statistical attributes
	m_AvgSpec.assign(numBands, 0.0f);
	m_StdSpec.assign(numBands, 0.0f);

	// Deserialize the average spectral values
	std::memcpy(&m_AvgSpec[0], &(stream[position]), numBands * FloatSize);
	position += numBands * FloatSize;

	// Deserialize the standard deviation values
	std::memcpy(&m_StdSpec[0], &(stream[position]), numBands * FloatSize);
	position += numBands * FloatSize;

	// Deserialize the perimeter
	std::memcpy(&m_Perimeter, &(stream[position]), UInt32Size);
	position += UInt32Size;

	// Deserialize the area
	std::memcpy(&m_Area, &(stream[position]), UInt32Size);
	position += UInt32Size;
}

}
}

#endif
