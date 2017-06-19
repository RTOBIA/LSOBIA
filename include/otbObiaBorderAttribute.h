#ifndef otbObiaBorderAttribute_h
#define otbObiaBorderAttribute_h

#include "otbObiaGenericAttribute.h"

/**
\file otbObiaBorderAttribute.h
\brief This file define the Mean Attribute class used to compute attribute of all features
*/
namespace otb
{
namespace obia
{


/** \class BorderAttribute
 *	\brief Class defining the MeanAttribute. It allows to compute the mean of exterior ring with a margin.
 */
template< class TInputImage >
class BorderAttribute : public GenericAttribute<TInputImage>
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self           = BorderAttribute;
    using Pointer        = itk::SmartPointer<Self>;
    using ConstPointer   = itk::SmartPointer< const Self>;
    using PixelType	     = typename TInputImage::InternalPixelType;
    using OGRLayerType 	 = ogr::Layer;
    using OGRFeatureType = ogr::Feature;

    BorderAttribute();
    BorderAttribute(double margin);
    virtual ~BorderAttribute();

    /**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
    virtual OGRFeatureType GenerateFeature(const OGRFeatureType inputFeature);

    /**Compute the attribut*/
    virtual double ComputeAttribut(std::vector<PixelType> samples);

    /**Set the margin*/
    void SetMargin(double margin){ m_Margin = margin;};

protected:

    /**Margin (default is 0)*/
    double m_Margin;

};

}
}
#include <otbObiaBorderAttribute.txx>
#endif
