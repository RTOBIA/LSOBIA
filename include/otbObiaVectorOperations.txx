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

bool
VectorOperations
::isFeatureValid(const OTBFeatureType feature)
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

std::vector<OGRLineString*>
VectorOperations
::SortLinesString(std::vector<OGRLineString*> unsortedGeoms)
{
	std::vector<OGRLineString*> sortedGeoms;
	std::vector<OGRLineString*> tmpGeoms = unsortedGeoms;
	bool isSorted = false;
	OGRPoint startingPoint, lastPoint;
	OGRLineString*  lastGeom = nullptr;
	if(unsortedGeoms.size() > 0)
	{
		//Initialize first element
		sortedGeoms.push_back(unsortedGeoms[0]);
		tmpGeoms.erase(tmpGeoms.begin(), tmpGeoms.begin() + 1);
		unsortedGeoms[0]->getPoint(0, &startingPoint);

		//Update last geom inserted
		lastGeom = sortedGeoms.back();
		lastGeom->getPoint(lastGeom->getNumPoints() - 1, &lastPoint);

		//std::cout << "First feature = " << unsortedGeoms[0]->exportToGML() << std::endl;
	}
	else
	{
		isSorted = true;
	}

	unsigned int cpt_while = 0;
	//While all not sorted
	while(!isSorted)
	{
		//Get last element of the sorted vector
		//Get the last point in the last element of the vector
		//Loop all features
		for(unsigned int fId = 0; fId < tmpGeoms.size(); ++fId)
		{
			//std::cout << "Fid = " << fId << "/" << tmpFeatures.size() << std::endl;
			//Check first coord
			OGRLineString* curGeom = tmpGeoms[fId];
			OGRPoint firstCurPoint, lastCurPoint;

			//First cur point
			curGeom->getPoint(0, &firstCurPoint);

			//Last cur point
			curGeom->getPoint(curGeom->getNumPoints() - 1, &lastCurPoint);

			//if points are equals
			if(lastPoint.Equals(&firstCurPoint))
			{
				//std::cout << "Point equals = " << lastPoint.exportToKML() << "/" << firstCurPoint.exportToGML() << std::endl;
				//Then add the edge after
				sortedGeoms.push_back(curGeom);

				//Remove from tmp vector
				tmpGeoms.erase(tmpGeoms.begin() + fId, tmpGeoms.begin() + fId + 1);

				//break the loop
				break;
			}
			else if(lastPoint.Equal(&lastCurPoint))
			{
				//std::cout << "Reverse Point equals = " << lastPoint.exportToKML() << "/" << firstCurPoint.exportToGML() << std::endl;
				//Then reverse the linestring
				curGeom->reversePoints();

				//Then add the edge after
				sortedGeoms.push_back(curGeom);

				//Remove from tmp vector
				tmpGeoms.erase(tmpGeoms.begin() + fId, tmpGeoms.begin() + fId + 1);
			}
		}

		//Update last geom
		lastGeom = sortedGeoms.back();
		lastGeom->getPoint(lastGeom->getNumPoints() - 1, &lastPoint);

		//DEBUG
		/*std::cout << "--------------- SORTED -----------" << std::endl;
		for(unsigned int l = 0; l < sortedGeoms.size(); ++l)
		{
			std::cout <<sortedGeoms[l]->exportToGML() <<std::endl;
		}
		std::cout << "---------------------------------" << std::endl;

		std::cout << "----------------REMANING-----------" << std::endl;
		for(unsigned int l = 0; l < tmpGeoms.size(); ++l)
		{
			std::cout <<tmpGeoms[l]->exportToGML() <<std::endl;
		}*/
		//if the last point of the last element inserted is equal to the first point, then we can close the external ring
		if(lastPoint.Equals(&startingPoint))
		{
			//std::cout << "Sorted and it remains : " << tmpGeoms.size() << std::endl;
			isSorted = true;
		}

		cpt_while++;
	}

	return sortedGeoms;
}

OGRPolygon*
VectorOperations
::CastToPolygon(OGRGeometry* geomRef)
{
	OGRwkbGeometryType geomType = geomRef->getGeometryType();
	if(geomType == wkbPolygon)
	{
		return (OGRPolygon*) geomRef;
	}
	else
	{
		std::cout << "Geometry " << geomRef->exportToKML() << " is not a polygon " << std::endl;
		return nullptr;
	}
}


//Create a field
void
VectorOperations
::CreateNewField(OTBLayerType& poLayer, std::string fieldName, OGRFieldType fieldType)
{
	OGRFieldDefn fieldDef(fieldName.c_str(), fieldType);
	poLayer.CreateField(fieldDef, true);

}

void
VectorOperations
::DisplayFeature(const OTBFeatureType feature, bool displayGeom)
{
	//Number of fields
	int fieldsCount = feature.GetDefn().GetFieldCount();
	for(unsigned fId = 0; fId < fieldsCount; ++fId)
	{
		std::cout << "---------------------------------------------------" << std::endl;
		OGRFieldDefn* fieldDef = feature.GetDefn().GetFieldDefn(fId);
		const char* fieldName = fieldDef->GetNameRef();
		OGRFieldType fieldType = fieldDef->GetType();

		std::cout << "Field " <<fieldDef->GetNameRef() << std::endl;
		std::cout <<"Field type = " << fieldDef->GetFieldTypeName(fieldType) << std::endl;
		//Get value
		switch(fieldType)
		{
			case OFTInteger64:
			{
				std::cout << "Value = " << feature.ogr().GetFieldAsDouble(fieldName) << std::endl;
				break;
			}
			case OFTInteger64List:
			{
				int nbValues = 0;
				const GIntBig* values = feature.ogr().GetFieldAsInteger64List(fieldName, &nbValues);
				for(int k = 0; k < nbValues; ++k)
				{
					std::cout << "Element " << k << " = " << values[k] << std::endl;;
				}
				break;
			}
			default:
			{
				break;
			}
		}
		std::cout << "---------------------------------------------------" << std::endl;
	}

	if(displayGeom)
	{
		std::cout << feature.ogr().GetGeometryRef()->exportToKML() << std::endl;
	}
}

bool
VectorOperations
::IsVerticalOrHorizontal(const OGRGeometry* geom)
{
	if(geom->getGeometryType() == wkbLineString)
	{
		//Get first point and last point
		OGRPoint firstPoint,lastPoint;
		OGRLineString* linestring = (OGRLineString*) geom;
		if(linestring->getNumPoints() == 2)
		{
			//If output linestring contains  2 point and both x or both y are equals, then the line is vertical or horizontal
			linestring->getPoint(0, &firstPoint);
			linestring->getPoint(1, &lastPoint);

			if(firstPoint.getX() == lastPoint.getX() ||
			   firstPoint.getY() == lastPoint.getY())
			{
				return true;
			}
		}
	}

	return false;



}

bool
VectorOperations
::IsRectangle(const OGRPolygon* polygon)
{
	//Get exterior ring
	const OGRLinearRing* extRing = polygon->getExteriorRing();

	//The 5th point correspond to the 1st point, so a rectangle is made with 5 points (the 4 corners, and the first point
	//in order to close the polygon)
	if(extRing->getNumPoints() != 5)
	{
		return false;
	}
	else
	{
		OGRPoint lastPoint;
		extRing->getPoint(0, &lastPoint);

		//Loop accross points
		double oldX = lastPoint.getX();
		double oldY = lastPoint.getY();

//		for(unsigned int i = 1; i < (extRing->getNumPoints() - 1); ++i)
//		{
//
//		}
		//Get bounding box
		OGRGeometry* boundary = polygon->Boundary();

		//If bounding box equal exterior ring, then polygon is a rectangle
		if(extRing->Equals(boundary))
		{
			std::cout << "isRectangle = " << polygon->exportToKML() << std::endl;
			return true;
		}
		else
		{
			return false;
		}

	}
}

std::vector<char>
VectorOperations
::SerializeLayer(const OTBLayerType layer)
{

}

VectorOperations::OTBLayerType
VectorOperations
::DeSerializeLayer(const std::vector<char>& serializedLayer)
{

}
} // end of namespace obia
} // end of namespace otb

#endif


