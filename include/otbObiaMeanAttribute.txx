#ifndef otbObiaMeanAttribute_txx
#define otbObiaMeanAttribute_txx

#include <otbObiaMeanAttribute.h>

namespace otb
{
namespace obia
{

template< class TInputImage >
MeanAttribute<TInputImage>
::MeanAttribute()
{
	this->m_FieldName 	 = "Mean";
	this->m_AttributName = "Mean Attribute";
	this->m_FieldType	 = OFTReal;
}

template< class TInputImage >
MeanAttribute<TInputImage>
::~MeanAttribute()
{

}

/**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
template< class TInputImage >
typename MeanAttribute<TInputImage>::OGRFeatureType
MeanAttribute<TInputImage>
::GenerateFeature(const OGRFeatureType inputFeature)
{
	//Create a new feature (clone of input feature because for the mean, we need the full geometry)
	OGRFeatureType requiredFeature(inputFeature);

	return requiredFeature;
}

/**Compute the attribut*/
template< class TInputImage >
double
MeanAttribute<TInputImage>
::ComputeAttribut(std::vector<PixelType> samples)
{
	double meanValue = 0.0;
	//Loop each samples
	for(unsigned int k = 0; k < samples.size(); k++)
	{
		PixelType sample = samples[k];
		meanValue += static_cast<double>(sample);
	}

	if(samples.size() != 0)
	{
		meanValue = meanValue/samples.size();
	}

	return meanValue;
}

}//end obia
}//end otb

#endif
