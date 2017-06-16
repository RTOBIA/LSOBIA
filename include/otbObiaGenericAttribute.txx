#ifndef otbObiaGenericAttribute_txx
#define otbObiaGenericAttribute_txx

#include <otbObiaGenericAttribute.h>

namespace otb
{
namespace obia
{

template< class TInputImage >
GenericAttribute<TInputImage>
::GenericAttribute()
{
	m_AttributName = "Generic Attribute";
	m_FieldName    = "GenericField";
	m_FieldType	   = OFTReal;
}

template< class TInputImage >
GenericAttribute<TInputImage>
::~GenericAttribute()
{
}


}//end obia
}//end otb

#endif
