#ifndef otbObiaMeanAttribut_h
#define otbObiaMeanAttribut_h

#include "otbObiaGenericAttribut.h"

namespace otb
{
namespace obia
{
template< class TInputImage >
class MeanAttribut : public GenericAttribut<TInputImage>
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self           = MeanAttribut;
    using Pointer        = itk::SmartPointer<Self>;
    using ConstPointer   = itk::SmartPointer< const Self>;
    using PixelType	     = typename TInputImage::InternalPixelType;
    using OGRLayerType 	 = ogr::Layer;
    using OGRFeatureType = ogr::Feature;

    /** Method for creation through the object factory. */
    MeanAttribut();
    virtual ~MeanAttribut();

    /**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
    virtual OGRFeatureType GenerateFeature(const OGRFeatureType inputFeature);

    /**Compute the attribut*/
    virtual double ComputeAttribut(std::vector<PixelType> samples);

};

}
}
#include "otbObiaMeanAttribut.txx"
#endif
