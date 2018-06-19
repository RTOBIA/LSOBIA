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
#ifndef otbObiaDouglasPeukerSimplify_txx
#define otbObiaDouglasPeukerSimplify_txx

#include "otbObiaDouglasPeukerSimplify.h"

namespace otb
{
namespace obia
{

DouglasPeukerFunc
::DouglasPeukerFunc() : GenericSimplifyFunc()
{
	m_tolerance = 1.0;
}

DouglasPeukerFunc
::~DouglasPeukerFunc()
{
	/**This class do not delete output geom.
	 * Pay attention with m_outputGeom to delete it after...*/
}


void
DouglasPeukerFunc
::SimplifyLine()
{
	/**Call GDAL Simplify for this geometry*/
	if(m_inputGeom->getGeometryType() == wkbMultiLineString)
	{
		//Force to linestring (clone the input. The memory is managed outside this class)
		//forceToLineString consumes the passed geometry, so by cloning we ensure to keep inputgeom
		m_outputGeom = OGRGeometryFactory::forceToLineString(m_inputGeom->clone(), false);

		//Care: forcetolinestring can still produce a multilinestring which cannot be simplified
		m_outputGeom = m_outputGeom->SimplifyPreserveTopology(this->m_tolerance);

		//NO MORE USEFULL SINCE WE CAN HANDLE INVALID GEOMETRIES
		//		//Check if the output line if vertical or horizontal
//		if(VectorOperations::IsVerticalOrHorizontal((OGRLineString*) m_outputGeom))
//		{
//			//Do not use this geometry, clone the entry one
//			delete m_outputGeom;
//			m_outputGeom = m_inputGeom->clone();
//		}
	}
	else if(m_inputGeom->getGeometryType() == wkbLineString ||
			m_inputGeom->getGeometryType() == wkbPoint		||
			m_inputGeom->getGeometryType() == wkbMultiPoint)
	{
		m_outputGeom = m_inputGeom->clone();
	}
	else
	{
		//TODO : theorically, we are not supposed to be there...We should have only multiline string and linestring
		//because of the intersection between geometries before going here
		std::cout << "Geometry not valid :: " << m_inputGeom->getGeometryName() << std::endl;
		std::cout << "Geometry not valid :: " << m_inputGeom->exportToGML() << std::endl;
		//std::cout << m_inputGeom->exportToKML() << std::endl;
		m_outputGeom = nullptr;
	}

	//std::cout <<" SIMPLIFIED ==> " << m_outputGeom->exportToKML() << std::endl;

}
}//end obia
}//end otb


#endif
