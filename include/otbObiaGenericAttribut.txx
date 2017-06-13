#ifndef otbObiaGenericAttribut_txx
#define otbObiaGenericAttribut_txx

#include "otbObiaGenericAttribut.h"

namespace otb
{
namespace obia
{

template< class TInputImage >
GenericAttribut<TInputImage>
::GenericAttribut()
{
	m_AttributName = "Generic Attribut";
	m_FieldName    = "GenericField";
	m_FieldType	   = OFTReal;
}

template< class TInputImage >
GenericAttribut<TInputImage>
::~GenericAttribut()
{
}


}//end obia
}//end otb

#endif
