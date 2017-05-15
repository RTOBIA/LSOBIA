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
	m_tolerance = 0.0;
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
	if(m_inputGeom->getGeometryType() == wkbMultiLineString)
	{
		//Force to linestring
		m_outputGeom = OGRGeometryFactory::forceToLineString(m_inputGeom, false);
		m_outputGeom = this->m_inputGeom->SimplifyPreserveTopology(this->m_tolerance);
		//exit(1);
	}else{
		m_outputGeom = m_inputGeom->clone();
	}

}
}//end obia
}//end otb


#endif
