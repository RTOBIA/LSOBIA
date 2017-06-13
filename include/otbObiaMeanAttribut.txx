#ifndef otbObiaMeanAttribut_txx
#define otbObiaMeanAttribut_txx

#include "otbObiaMeanAttribut.h"

namespace otb
{
namespace obia
{

template< class TInputImage >
MeanAttribut<TInputImage>
::MeanAttribut()
{
	this->m_FieldName 	 = "Mean";
	this->m_AttributName = "Mean Attribut";
	this->m_FieldType	 = OFTReal;
}

template< class TInputImage >
MeanAttribut<TInputImage>
::~MeanAttribut()
{

}

/**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
template< class TInputImage >
typename MeanAttribut<TInputImage>::OGRFeatureType
MeanAttribut<TInputImage>
::GenerateFeature(const OGRFeatureType inputFeature)
{
	//Create a new feature (clone of input feature because for the mean, we need the full geometry)
	OGRFeatureType requiredFeature(inputFeature);

	return requiredFeature;
}

/**Compute the attribut*/
template< class TInputImage >
double
MeanAttribut<TInputImage>
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
