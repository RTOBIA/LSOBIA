#ifndef otbObiaVectorOperations_txx
#define otbObiaVectorOperations_txx

#include "otbObiaVectorOperations.h"
#include <stdio.h>
#include <string.h>
namespace otb
{
namespace obia
{
template<typename T>
VectorOperations::OTBFeatureType
VectorOperations::GetFeatureAt(OGRDataSourceType* ogrDS,
   		  	  	  	  	   	   const std::string layerName,
							   const std::string fieldName,
							   const T fieldValue)

{
	//Get the layer
	OTBLayerType layer = ogrDS->GetLayer(layerName);
	//Loop accross all feature
	//unsigned int nbFeatures = layer.GetFeatureCount(true);
	unsigned int nbFeatures = layer.GetFeatureCount(true);
	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		//Current feature
		OTBFeatureType curFeature = layer.GetFeature(fId);

		//Maybe usefull another time
		OGRFieldType fieldType = curFeature.GetFieldDefn(fieldName).GetType();

		//Extract field
		unsigned int curValue = curFeature.ogr().GetFieldAsInteger(fieldName.c_str());

		//Check field
		if(curValue == fieldValue)
		{
			//If the value we are looking for is ok, return the feature
			return curFeature;
		}

	}


	return nullptr;

}


template<typename T>
std::vector<VectorOperations::OTBFeatureType>
VectorOperations
::GetAllFeaturesAt(OGRDataSourceType* ogrDS,
   		  	  	   const std::string layerName,
				   const std::string fieldName,
				   const T fieldValue)

{

}

std::vector<VectorOperations::OTBFeatureType>
VectorOperations
::GetAdjacentFeatures(OGRDataSourceType* ogrDS,
	  	  	   	   	  const std::string layerName,
					  const OTBFeatureType feature)
{
	//Initialise vector
	std::vector<OTBFeatureType> adjFeatures;

	//Extract coordinates of all adjacents feature
	int nbAdjacents = 0;
	const long long int* adjCoordinates = feature.ogr().GetFieldAsInteger64List(adjStartingCoordsFieldName.c_str(),
																		  &nbAdjacents);

	//Loop accross all adjacent coordiantes
	for(unsigned int fId = 0; fId < nbAdjacents; ++fId)
	{
		//Get feature
		OTBFeatureType adjFeature = GetFeatureAt<unsigned int>(ogrDS, layerName, startingCoordsFieldName, adjCoordinates[fId]);
		adjFeatures.push_back(adjFeature);
	}

	return adjFeatures;
}

bool VectorOperations::isFeatureValid(const OTBFeatureType feature)
{
	//Get field value
	std::string fieldValue = feature.ogr().GetFieldAsString(validPolygonFieldName.c_str());
	if(strcmp(fieldValue.c_str(), "true") == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<char> VectorOperations::SerializeLayer(const OTBLayerType layer)
{

}

VectorOperations::OTBLayerType VectorOperations::DeSerializeLayer(const std::vector<char>& serializedLayer)
{

}
} // end of namespace obia
} // end of namespace otb

#endif

