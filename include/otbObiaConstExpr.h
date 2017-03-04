#ifndef otbObiaConstExpr_h
#define otbObiaConstExpr_h
#include <cstdint>

namespace otb
{

namespace obia
{

/** Some convenient alias for the whole obia chain */
using IdType = uint32_t;
using CoordValueType = uint32_t;

constexpr uint64_t BoolSize       = sizeof(bool);
constexpr uint64_t CharSize       = sizeof(char);
constexpr uint64_t IntSize        = sizeof(int);
constexpr uint64_t UInt8Size      = sizeof(uint8_t);
constexpr uint32_t UInt32Size     = sizeof(uint32_t);
constexpr uint64_t UInt64Size     = sizeof(uint64_t);
constexpr uint64_t FloatSize      = sizeof(float);
constexpr uint64_t IdSize         = sizeof(IdType);
constexpr uint64_t CoordValueSize = sizeof(CoordValueType);

} // end of namespace obia

} // end of namespace otb

#endif