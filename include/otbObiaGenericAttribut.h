#ifndef otbObiaGenericAttribut_h
#define otbObiaGenericAttribut_h

#include <otbOGRDataSourceWrapper.h>

namespace otb
{
namespace obia
{
template< class TInputImage >
class GenericAttribut
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self           = GenericAttribut;
    using Pointer        = itk::SmartPointer<Self>;
    using ConstPointer   = itk::SmartPointer< const Self>;
    using PixelType	     = typename TInputImage::InternalPixelType;
    using OGRLayerType 	 = ogr::Layer;
    using OGRFeatureType = ogr::Feature;

    /** Method for creation through the object factory. */
    GenericAttribut();
    virtual ~GenericAttribut();

    /**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
    virtual OGRFeatureType GenerateFeature(const OGRFeatureType inputFeature) = 0;

    /**Compute the attribut*/
    virtual double ComputeAttribut(std::vector<PixelType> samples) = 0;

    //Get attribut name
    std::string GetAttributName(){return m_AttributName;};

    //Get field Name
    std::string GetFieldName(){return m_FieldName;};

    //Get field type
    OGRFieldType GetFieldType(){return m_FieldType;};

protected:

    /**Field name*/
    std::string m_FieldName;

    //Attribut name
    std::string m_AttributName;

    //Field type
    OGRFieldType m_FieldType;
};

}
}
#include "otbObiaGenericAttribut.txx"
#endif

