#ifndef otbObiaVectorOperations_txx
#define otbObiaVectorOperations_txx

#include "otbObiaVectorOperations.h"
#include <stdio.h>
#include <string.h>
#include "otbMPIConfig.h"
#include "otbObiaConstExpr.h"

namespace otb
{
namespace obia
{
template<typename T>
VectorOperations::OGRFeatureType
VectorOperations::GetFeatureAt(OGRDataSourceType* ogrDS,
   		  	  	  	  	   	   const std::string layerName,
							   const std::string fieldName,
							   const T fieldValue)

{
	//Get the layer
	OGRLayerType layer = ogrDS->GetLayer(layerName);
	//Loop accross all feature
	//unsigned int nbFeatures = layer.GetFeatureCount(true);
	unsigned int nbFeatures = layer.GetFeatureCount(true);
	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		//Current feature
		OGRFeatureType curFeature = layer.GetFeature(fId);

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
std::vector<VectorOperations::OGRFeatureType>
VectorOperations
::GetAllFeaturesAt(OGRDataSourceType* ogrDS,
   		  	  	   const std::string layerName,
				   const std::string fieldName,
				   const T fieldValue)

{

}

std::vector<VectorOperations::OGRFeatureType>
VectorOperations
::GetAdjacentFeatures(OGRDataSourceType* ogrDS,
	  	  	   	   	  const std::string layerName,
					  const OGRFeatureType feature)
{
	//Initialise vector
	std::vector<OGRFeatureType> adjFeatures;

	//Extract coordinates of all adjacents feature
	int nbAdjacents = 0;
	const long long int* adjCoordinates = feature.ogr().GetFieldAsInteger64List(adjStartingCoordsFieldName.c_str(),
																				&nbAdjacents);

	//std::cout << feature.GetGeometry()->exportToGML() << std::endl;

	//Loop accross all adjacent coordiantes
	for(unsigned int fId = 0; fId < nbAdjacents; ++fId)
	{
		//std::cout << "Looking for feature at " << adjCoordinates[fId] << std::endl;
		//Get feature
		OGRFeatureType adjFeature = VectorOperations::GetFeatureAt<unsigned int>(ogrDS, layerName, startingCoordsFieldName, adjCoordinates[fId]);
		adjFeatures.push_back(adjFeature);
		//std::cout << "Adjacent feature " << fId << " = \n" << adjFeature.GetGeometry()->exportToGML() << std::endl;
	}

	//Sort feature according to their coordinates to get the first feature with the lowest coordinate
	std::sort(adjFeatures.begin(), adjFeatures.end(), [](const OGRFeatureType a,
													   const OGRFeatureType b)->bool{
			double aCoord = a.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
			double bCoord = b.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
			return aCoord < bCoord;
	});


	return adjFeatures;
}

bool
VectorOperations
::isFeatureValid(const OGRFeatureType feature)
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
		//std::cout << "CPT WHILE = " << cpt_while << std::endl;
		//Get last element of the sorted vector
		//Get the last point in the last element of the vector
		//Loop all features
		for(unsigned int fId = 0; fId < tmpGeoms.size(); ++fId)
		{
			//std::cout << "Fid = " << fId << "/" << tmpFeatures.size() << std::endl;
			//Check first coord
			OGRLineString* curGeom = tmpGeoms[fId];
			OGRPoint firstCurPoint, lastCurPoint;

			//std::cout << "GEOMETRY " << fId << " = " << curGeom->exportToGML() << std::endl;
			//First cur point
			curGeom->getPoint(0, &firstCurPoint);

			//Last cur point
			curGeom->getPoint(curGeom->getNumPoints() - 1, &lastCurPoint);

			//if points are equals
			//std::cout << "Compare " << lastPoint.exportToGML() << " and " << firstCurPoint.exportToGML() << std::endl;
			//std::cout << "Compare " << lastPoint.exportToGML() << " and " << lastCurPoint.exportToGML() << std::endl;
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
			else if(lastPoint.Equals(&lastCurPoint))
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

		//isSorted = true;

		//Update last geom
		lastGeom = sortedGeoms.back();
		lastGeom->getPoint(lastGeom->getNumPoints() - 1, &lastPoint);

		//if the last point of the last element inserted is equal to the first point, then we can close the external ring
		if(lastPoint.Equals(&startingPoint))
		{
			//std::cout << "Sorted and it remains : " << tmpGeoms.size() << std::endl;
			isSorted = true;
		}

		cpt_while++;
		if(cpt_while > unsortedGeoms.size() + 1)
		{
			std::cout << "Impossible de trier pour " << MPIConfig::Instance()->GetMyRank() << std::endl;
			std::cout << "============= UNSORTED ===============" << std::endl;
			for(int i = 0; i < unsortedGeoms.size() ; ++i)
			{
				std::cout << "  <gml:featureMember> " << std::endl;
				std::cout << "\t\t <ogr:cleaned fid=\"cleaned.0\"> " << std::endl;
				std::cout << unsortedGeoms[i]->exportToGML() << std::endl;

				std::cout << " </ogr:cleaned> "  << std::endl;
				std::cout << "\t\t   </gml:featureMember> " << std::endl;

			}

			std::cout << "============= SORTED ===============" << std::endl;
			for(int i = 0; i < sortedGeoms.size(); ++i)
			{
				std::cout << sortedGeoms[i]->exportToGML() << std::endl;
			}
			sortedGeoms.clear();
			return sortedGeoms;
		}
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
::CreateNewField(OGRLayerType& poLayer, std::string fieldName, OGRFieldType fieldType)
{
	OGRFieldDefn fieldDef(fieldName.c_str(), fieldType);
	poLayer.CreateField(fieldDef, true);

}

void
VectorOperations
::DisplayFeature(const OGRFeatureType feature, bool displayGeom)
{
	std::cout << "---------------------------------------------------" << std::endl;
	//Number of fields
	int fieldsCount = feature.GetDefn().GetFieldCount();
	for(unsigned fId = 0; fId < fieldsCount; ++fId)
	{
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
	}

	if(displayGeom)
	{
		std::cout << "Feature's geometry = " << std::endl;
		std::cout << feature.ogr().GetGeometryRef()->exportToGML() << std::endl;
	}

	std::cout << "---------------------------------------------------" << std::endl;
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
::IsParallelogram(const OGRPolygon* polygon)
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
		//Compute 4 length
		std::vector<double> ringLengths;
		for(unsigned int i = 0; i < (extRing->getNumPoints() - 1); ++i)
		{
			OGRPoint firstPoint,lastPoint;
			extRing->getPoint(i, &firstPoint);
			extRing->getPoint(i + 1, &lastPoint);

			//Compute distance
			double d = firstPoint.Disjoint(&lastPoint);
			ringLengths.push_back(d);
		}

		//Compare lenght 0 and 2, and lenght 1 and 3
		if(fabs(ringLengths[0] - ringLengths[2]) < std::numeric_limits<double>::epsilon() &&
		   fabs(ringLengths[1] - ringLengths[3]) < std::numeric_limits<double>::epsilon() )

		{
			return true;
		}
		else
		{
			return false;
		}
//		//Get bounding box
//		OGRGeometry* boundary = polygon->Boundary();
//
//		//If bounding box equal exterior ring, then polygon is a rectangle
//		if(extRing->Equals(boundary))
//		{
//			return true;
//		}
//		else
//		{
//			return false;
//		}

	}
}


bool
VectorOperations
::IsValid(const OGRGeometry* geom)
{
	OGRwkbGeometryType geomType = geom->getGeometryType();
	switch(geomType)
	{
		case wkbPolygon:
		{
			const OGRPolygon* poly = static_cast<const OGRPolygon*> (geom);
			return IsValid(poly->getExteriorRing());
		}
		default:
		{
			return geom->IsValid();
		}
	}
}

std::vector<OGRGeometry*>
VectorOperations
::IntersectGeoms(OGRGeometry* geomRef, OGRGeometry* geomAdj)
{
	std::vector<OGRGeometry*> edges;

	OGRGeometry* intersectedLine = nullptr;

	if(geomRef->Touches(geomAdj))
	{
		intersectedLine = geomRef->Intersection(geomAdj);
		OGRwkbGeometryType geomType = intersectedLine->getGeometryType();

		switch(geomType)
		{
			case wkbGeometryCollection:
			case wkbMultiLineString:
			{
				OGRGeometryCollection* geomCollec = (OGRGeometryCollection*) intersectedLine;
				unsigned int nbGeoms = geomCollec->getNumGeometries();
				for(unsigned int geomId = 0; geomId < nbGeoms; ++geomId)
				{
					edges.push_back(geomCollec->getGeometryRef(geomId));
				}
				break;
			}
			default:
				edges.push_back(intersectedLine);
				break;
		}
	}//End if touch
	return edges;
}
GDALDriver*
VectorOperations
::InitializeGDALDriver(std::string driverName)
{
	GDALDriver *poDriver;
	GDALAllRegister();
	poDriver = (GDALDriver*) GDALGetDriverByName(driverName.c_str());
	if( poDriver == NULL )
	{
		printf( "%s driver not available\n", driverName.c_str());
		exit( 1 );
	}

	return poDriver;

}

GDALDataset*
VectorOperations
::CloneDataset(GDALDataset* sourceDs)
{
	std::cout << "-----------------Clone Dataset ---------------------" << std::endl;
	//Intialize a driver
	GDALDriver* poDriver = InitializeGDALDriver(sourceDs->GetDriverName());

	//Initialize dataset
	GDALDataset* poDs = poDriver->CreateCopy("", sourceDs, true, nullptr, nullptr, nullptr);

	//Return created
	return poDs;


}

void
VectorOperations
::WriteOGRDataSource(OGRDataSourceType* ogrDS, std::string filename, int layerId)
{
	//Get output
	//OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );
	GDALDriver*  poDriverGml = VectorOperations::InitializeGDALDriver("GML");
	GDALDataset* polyGDs = poDriverGml->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );

	if(layerId == -1)
	{
		//Write all layers
		for(unsigned int layerId = 0; layerId < ogrDS->ogr().GetLayerCount(); ++layerId)
		{
			std::cout << "--------------------------------------------------------" << std::endl;
			std::cout << "Writing layer " << ogrDS->GetLayer(layerId).ogr().GetName() << std::endl;
			std::cout << "Number feature = " << ogrDS->GetLayer(layerId).ogr().GetFeatureCount(true) << std::endl;
			std::cout << "--------------------------------------------------------" << std::endl;
			OGRLayerType curLayer = ogrDS->GetLayer(layerId);
			polyGDs->CopyLayer(&(curLayer.ogr()), curLayer.GetName().c_str());
		}
	}
	else
	{
		otb::ogr::Layer curLayer = ogrDS->GetLayer(layerId);
		polyGDs->CopyLayer(&(curLayer.ogr()), curLayer.GetName().c_str());
	}


	//polyGDs->CopyLayer(&(layer.ogr()), layer.GetName().c_str(), nullptr);
	//Clean memory
	//Maybe see if we need to delete driver...
	GDALClose(polyGDs);
}

void
VectorOperations
::WriteOGRDataSource(OGRDataSourceType* ogrDS, std::string filename, std::string layerName)
{
	if(layerName.compare("") == 0)
	{
		WriteOGRDataSource(ogrDS, filename, -1);
	}
	else
	{
		//Loop layer
		for(unsigned int layerId = 0; layerId < ogrDS->ogr().GetLayerCount(); ++layerId)
		{
			std::string currentLayerName = ogrDS->GetLayer(layerId).GetName();
			if(currentLayerName.compare(layerName))
			{
				WriteOGRDataSource(ogrDS, filename, layerId);
				return;
			}
		}

	}
}
//
OGRGeometry* VectorOperations
::ConvertToMultiPolygons(OGRPolygon* polygon)
{
	std::cout << "CONVERT TO MULTIPOLYGONS" << std::endl;
	std::cout << "=========================================================" << std::endl;
	std::cout << polygon->exportToGML() << std::endl;
	std::cout << "=========================================================" << std::endl;
  	//Convert exterior ring to multilinestring
	std::vector<OGRLineString*> polygonEdges = ConvertPolygonToLinearString(polygon);

	//Loop through each edges, and compute intersections
	for(int edgeId = 0; edgeId < polygonEdges.size(); ++edgeId)
	{
		//std::cout << "-----------------------------------------" << std::endl;
		//Current linestring
		OGRLineString* curLinestring = polygonEdges[edgeId];
		//std::cout << "CUR LINE = " << curLinestring->exportToGML() << std::endl;

		//Intersect with all others
		for(int l = 0; l < polygonEdges.size(); ++l)
		{
			OGRLineString* otherLinestring = polygonEdges[l];
			if(edgeId != l)
			{
//				std::cout << "-----------------------------------------" << std::endl;
//				std::cout << "OTHER LINE = " << otherLinestring->exportToGML() << std::endl;

				//Compute intersection (must be a point because each linestring is composed by 2 points)
				if(curLinestring->Intersects(otherLinestring))
				{
					OGRGeometry* intersectedPoint = curLinestring->Intersection(otherLinestring);
					//std::cout << "POINT = " << intersectedPoint->exportToGML() << std::endl;

					//if it is not a point, it means that the lines are overlapping, we do nothing
					if(intersectedPoint->getGeometryType() == wkbPoint)
					{
						InsertPointInRightOrder(curLinestring, (OGRPoint*) intersectedPoint);
						InsertPointInRightOrder(otherLinestring, (OGRPoint*) intersectedPoint);
						//std::cout << "-----------------------------------------" << std::endl;
						//std::cout << "NEW OTHER LINE STRING = "<< otherLinestring->exportToGML() << std::endl;

					}
				}

				//std::cout << "-----------------------------------------" << std::endl;

			}
		}
//		std::cout << "NEW LINE = " << curLinestring->exportToGML() << std::endl;
//		std::cout << "-----------------------------------------" << std::endl;

	}

	//Create a multilinestring
	OGRMultiLineString* multilinestring = new OGRMultiLineString();
	for(unsigned int i = 0; i < polygonEdges.size(); ++i)
	{
		OGRLineString* curLinestring = polygonEdges[i];
		for(int j = 0; j < (curLinestring->getNumPoints() - 1); ++j)
		{
			OGRLineString* tmp = new OGRLineString();
			OGRPoint firstPoint, lastPoint;
			curLinestring->getPoint(j, &firstPoint);
			curLinestring->getPoint(j + 1, &lastPoint);
			tmp->addPoint(&firstPoint);
			tmp->addPoint(&lastPoint);

			multilinestring->addGeometryDirectly(tmp);
		}
		//Delete the old geom
		delete curLinestring;
	}

	//Polygonize
	OGRGeometry*  multipolygons = multilinestring->Polygonize();

	if(multipolygons == nullptr)
	{
		std::cout << "FAILED TO POLYGONIZE : " << multilinestring->exportToGML() << std::endl;
	}

	return multipolygons;
}
void VectorOperations
::InsertPointInRightOrder(OGRLineString* linestring, OGRPoint* point)
{
	//std::cout << "POINT = " << point->exportToGML() << std::endl;

	//If the point is equal to end or start point, then do nothing
	OGRPoint startPoint;
	OGRPoint lastPoint;

	//Compute distance from start point to point
	linestring->StartPoint(&startPoint);
	linestring->EndPoint(&lastPoint);

	//First: we check if this point is a first or last point of the linestring
	//meaning that the intersection is between adjacent linestring
	if(point->Equals(&startPoint) || point->Equals(&lastPoint))
	{
		return;
	}

	//If not, and the linestring contains only 2 points, we can add safely this point between
	if(linestring->getNumPoints() < 3)
	{
		//std::cout << "Inserting " << point->exportToGML() << " before " << lastPoint.exportToGML() << std::endl;
 		linestring->setPoint(1, point);
		linestring->setPoint(2, &lastPoint);
		return;
	}

	double distance = startPoint.Distance(point);

	//std::cout << "LINE STRING INSERTING POINT : " << linestring->exportToGML() <<std::endl;
	//Loop each point (except first)
	for(int i = 1; i < linestring->getNumPoints(); ++i)
	{
		OGRPoint curPoint;
		linestring->getPoint(i, &curPoint);

		if(!curPoint.Equals(point))
		{
			//Compute distance from start point
			double curDistance = startPoint.Distance(&curPoint);

			//If distance < curDistance, then we insert point before this one
			if(distance < curDistance)
			{
				//Set new point to i-position
				linestring->setPoint(i, point);

				//Shift all remaining points to the right
				for(unsigned int j = i + 1; j < linestring->getNumPoints(); ++j)
				{
					//Old point
					OGRPoint prevPoint;
					linestring->getPoint(j, &prevPoint);

					//Set new point
					linestring->setPoint(j, &curPoint);

					//Set new cur point
					curPoint = prevPoint;

				}

				//Set last point
				linestring->addPoint(&curPoint);

				//std::cout << "Insert at " << i << std::endl;

				return;
			}

		}
		else
		{
			//The point already exist, we just return
			return;
		}
	}
}
std::vector<OGRLineString*> VectorOperations
::ConvertPolygonToLinearString(OGRPolygon* polygon)
{
	OGRLinearRing* extRing = polygon->getExteriorRing();

	std::vector<OGRLineString*> polygonEdges;
	//Loop through each points and create successive linestring into a vector
	for(int pid = 0; pid < (extRing->getNumPoints() - 1); ++pid)
	{
		OGRPoint curPoint;
		OGRPoint nextPoint;

		//Get point and next point
		extRing->getPoint(pid, &curPoint);
		extRing->getPoint(pid + 1, &nextPoint);

		//Create linestring
		OGRLineString* linestring = new OGRLineString();
		linestring->addPoint(curPoint.getX(), curPoint.getY());
		linestring->addPoint(nextPoint.getX(), nextPoint.getY());

		//Add to vector
		polygonEdges.push_back(linestring);
	}

	return polygonEdges;
}
std::vector<char>
VectorOperations
::SerializeLayer(const OGRLayerType layer)
{

}

VectorOperations::OGRLayerType
VectorOperations
::DeSerializeLayer(const std::vector<char>& serializedLayer)
{

}
} // end of namespace obia
} // end of namespace otb

#endif


