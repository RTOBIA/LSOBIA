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

//Layer name
const std::string originalLayerName				= "original";
const std::string cleanedLayerName				= "cleaned";
const std::string nodataLayerName				= "nodata";
//Field for vector
const std::string labelFieldName		  	 	= "Label";
const std::string startingCoordsFieldName 		= "StartingCoords";
const std::string adjStartingCoordsFieldName 	= "AdjStartingCoords";
const std::string validPolygonFieldName		  	= "isValid";
const std::string reconstructedLayerName		= "Reconstructed";
const std::string unvalidLayerName				= "Unvalid";
const std::string attributesLayerName			= "AttributesLayer";

//Field for edge
const std::string polygonEdgeFieldName			= "PolygonCoords";

} // end of namespace obia

} // end of namespace otb

#endif
