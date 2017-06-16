#ifndef otbObiaGenericAttribute_h
#define otbObiaGenericAttribute_h

#include <otbOGRDataSourceWrapper.h>

/**
\file otbObiaMeanAttribute.h
\brief This file define the Generic Attribute class used to define all wanted attributes
*/
namespace otb
{
namespace obia
{
template< class TInputImage >
class GenericAttribute
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self           = GenericAttribute;
    using Pointer        = itk::SmartPointer<Self>;
    using ConstPointer   = itk::SmartPointer< const Self>;
    using PixelType	     = typename TInputImage::InternalPixelType;
    using OGRLayerType 	 = ogr::Layer;
    using OGRFeatureType = ogr::Feature;

    /**\brief Constructor*/
    GenericAttribute();

    /**\brief Destructor*/
    virtual ~GenericAttribute();

    /**\brief Generate the required feature (for example instead of full polygon, we could only need exterior ring)
     * This method is abstract, each attribute must implement it
     * \param : Input feature
     * \return : Required feature on which the attribute will be computed */
    virtual OGRFeatureType GenerateFeature(const OGRFeatureType inputFeature) = 0;

    /**\brief Compute attribute according to the samples given.
     * This method is abstract, each attribute must implement it
     * \param : Pixels from the ROI
     * \return : Computed attribute*/
    virtual double ComputeAttribut(std::vector<PixelType> samples) = 0;

    /**\brief Getter of attribute name
     * \return : Attribute name*/
    std::string GetAttributName(){return m_AttributName;};

    /**\brief Getter of field name
     * \return : Field name*/
    std::string GetFieldName(){return m_FieldName;};

    /**\brief Getter of field type
     * \return : Field type*/
    OGRFieldType GetFieldType(){return m_FieldType;};

protected:

    /**\brief Field name used in feature*/
    std::string m_FieldName;

    /**\brief Attribut name, used log*/
    std::string m_AttributName;

    /**\brief Field type used in feature to define the field*/
    OGRFieldType m_FieldType;
};

}
}
#include <otbObiaGenericAttribute.txx>
#endif

