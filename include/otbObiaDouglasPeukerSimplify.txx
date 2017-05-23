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
	 * Pay attentio nwith m_outputGeom to delete it after...*/
}


void
DouglasPeukerFunc
::SimplifyLine()
{
	/**Call GDAL Simplify for this geometry*/
	/**TODO*/
	//Check GDAL version, because the simplify preserve topology is undefined...
	std::cout <<" SIMPLIFY ==> " << m_inputGeom->exportToKML() << std::endl;
	if(m_inputGeom->getGeometryType() == wkbMultiLineString)
	{
		//Force to linestring (clone the input. The memory is managed outside this class)
		//forceToLineString consumes the passed geometry, so by cloning we nesure to keep inputgeom
		m_outputGeom = OGRGeometryFactory::forceToLineString(m_inputGeom->clone(), false);

		//Care: forcetolinestring can still produce a multilinestring which cannot be simplified
		//TODO: maybe add so robustness to that

		std::cout << "Line string : " << m_outputGeom->exportToKML() << std::endl;
		m_outputGeom = m_outputGeom->SimplifyPreserveTopology(this->m_tolerance);

		//Check if the output line if vertical or horizontal
		if(VectorOperations::IsVerticalOrHorizontal((OGRLineString*) m_outputGeom))
		{
			//Do not use this geometry, clone the entry one
			delete m_outputGeom;
			m_outputGeom = m_inputGeom->clone();
		}

		//std::cout << "AFTER OUTPUT GEOM = " << m_outputGeom->exportToGML() << std::endl;

		//exit(1);
	}
	else if(m_inputGeom->getGeometryType() == wkbLineString || m_inputGeom->getGeometryType() == wkbPoint)
	{
		m_outputGeom = m_inputGeom->clone();
	}
	else
	{
		//TODO : theorically, we are not supposed to be there...We should have only multiline string and linestring
		//because of the intersection between geometries before going here
		std::cout << "Geometry not valid :: " << m_inputGeom->getGeometryName() << std::endl;
		std::cout << m_inputGeom->exportToKML() << std::endl;
		m_outputGeom = nullptr;

		exit(EXIT_FAILURE);
	}

	std::cout <<" SIMPLIFIED ==> " << m_outputGeom->exportToKML() << std::endl;

}
}//end obia
}//end otb


#endif
