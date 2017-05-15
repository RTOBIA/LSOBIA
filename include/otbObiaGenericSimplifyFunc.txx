#ifndef otbObiaGenericSimplifyFunc_txx
#define otbObiaGenericSimplifyFunc_txx

#include "otbObiaGenericSimplifyFunc.h"

namespace otb
{
namespace obia
{

GenericSimplifyFunc
::GenericSimplifyFunc()
{
	m_inputGeom  = nullptr;
	m_outputGeom = nullptr;
}

GenericSimplifyFunc
::~GenericSimplifyFunc()
{
	/**This class do not delete output geom.
	 * Pay attention nwith m_outputGeom to delete it after...*/
}


}//end obia
}//end otb

#endif
