#ifndef otbObiaSimplifyVectorFilter_txx
#define otbObiaSimplifyVectorFilter_txx

#include "otbObiaSimplifyVectorFilter.h"
#include "otbObiaVectorOperations.h"

namespace otb
{
namespace obia
{

template <class TSimplifyFunc>
SimplifyVectorFilter<TSimplifyFunc>
::SimplifyVectorFilter() : m_FieldName("DN"), m_EdgeLayer(nullptr, true), m_NodataLayer(nullptr, true),
 m_BBFeature(nullptr), m_simplifyFunc(nullptr), m_CurrentCoords(0), m_EdgeFeatureDefn(nullptr)
{
   this->SetNumberOfRequiredInputs(1);
   this->SetNumberOfRequiredOutputs(1);

   GDALAllRegister();

   this->ProcessObject::SetNthOutput(0, this->MakeOutput(0) );
}

template <class TSimplifyFunc>
SimplifyVectorFilter<TSimplifyFunc>
::~SimplifyVectorFilter()
{
	//TODO
	//Free all maps
}
template <class TSimplifyFunc>
typename SimplifyVectorFilter<TSimplifyFunc>::DataObjectPointer
SimplifyVectorFilter<TSimplifyFunc>
::MakeOutput(DataObjectPointerArraySizeType itkNotUsed(idx))
{
  return static_cast< DataObjectPointer >(OGRDataSourceType::New().GetPointer());
}

template <class TSimplifyFunc>
const typename SimplifyVectorFilter<TSimplifyFunc>::OGRDataSourceType *
SimplifyVectorFilter<TSimplifyFunc>
::GetOutput()
{
  return static_cast< const OGRDataSourceType *>(
              this->ProcessObject::GetOutput(0));
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::SetInput(const OGRDataSourceType *input)
{
  this->Superclass::SetNthInput(0, const_cast<OGRDataSourceType *>(input));
}

template <class TSimplifyFunc>
const typename SimplifyVectorFilter<TSimplifyFunc>
::OGRDataSourceType *
SimplifyVectorFilter<TSimplifyFunc>
::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
    {
    return ITK_NULLPTR;
    }

  return static_cast<const OGRDataSourceType *>(this->Superclass::GetInput(0));
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::GenerateInputRequestedRegion(void)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename OGRDataSourceType::Pointer input  =
    const_cast<OGRDataSourceType *> (this->GetInput());

  if ( !input )
    {
    return;
    }
  // The input is necessarily the largest possible region.
  input->SetRequestedRegionToLargestPossibleRegion();

}


template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::GenerateData(void)
{
	std::cout << " --------------- SIMPLIFY VECTOR FILTER FOR "
			  << MPIConfig::Instance()->GetMyRank() <<"----------" << std::endl;
	//Get the input
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );

	//Initialize output
	InitializeOutputDS();

	//Loop accross all geometry
	std::cout << "Generate Data for SimplifyVectorFilter" << std::endl;

	//Get number of features
	OGRLayerType layer = ogrDS->GetLayer(m_LayerName);

	//Get the nodata layer
	m_NodataLayer = ogrDS->GetLayer(nodataLayerName);

	//Compute nb of features
	unsigned int nbFeatures = layer.GetFeatureCount(true);

	std::cout << "Nb features = " << nbFeatures << std::endl;

	//Creater bounding feature (bouding all the features)
	m_BBFeature = CreateBoundingBoxFeature();

	//add all features to a layer, and write it
	m_EdgeLayer = ogrDS->CreateLayer("Edges"  , ITK_NULLPTR, wkbUnknown);
	VectorOperations::CreateNewField(m_EdgeLayer, polygonEdgeFieldName, OFTInteger64List);

	//Set definition used for next feature
	m_EdgeFeatureDefn = &(m_EdgeLayer.GetLayerDefn());

	//Maybe not usefull
	OGRLayerType insidePolygons	= ogrDS->CreateLayer("Inside"  , ITK_NULLPTR, wkbUnknown);
	//VectorOperations::CreateNewField(insidePolygons, "Inside", OFTInteger64);
	//VectorOperations::CreateNewField(insidePolygons, "Englobing", OFTInteger64);

	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		std::cout << "Feature id " << fId << "/" << nbFeatures - 1 << std::endl;
		//Current feature
		OGRFeatureType curFeature = layer.GetFeature(fId);

		//Compute intersection with BB
		m_CurrentCoords = curFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
		m_bbEdges[m_CurrentCoords] = IntersectWithBoundaries(curFeature);

		if(m_CurrentCoords == 517731 && MPIConfig::Instance()->GetMyRank() == 0)
		{
			std::cout << "PROBLEM WITH " << curFeature.GetGeometry()->exportToGML() << std::endl;
		}

		//Initialize is englobed map value
		if(m_IsEnglobedMap.find(m_CurrentCoords) == m_IsEnglobedMap.end())
		{
			m_IsEnglobedMap[m_CurrentCoords] = false;
		}


		//Get all adjacent feature to the current one
		std::vector<OGRFeatureType> adjFeatures = VectorOperations::GetAdjacentFeatures(ogrDS,
																						m_LayerName, curFeature);

		//Once we got all adjacent features, we have to do the intersection between each adj and feature
		ConvertToEdges(curFeature, adjFeatures, insidePolygons);

		//Set false to is Valid field to not use it again
		curFeature.ogr().SetField(validPolygonFieldName.c_str(), "false");

		//Update feature to modify field value
		layer.SetFeature(curFeature);

//		//Create all edge feature to add into layer
//		for(unsigned int fId = 0; fId < edges.size(); ++fId)
//		{
//			//edgesLayer.SetFeature(*edges[fId]);
//			edgesLayer.CreateFeature(*edges[fId]);
//		}
	}

	//Get output pointer
	OGRDataSourceType::Pointer outputOGRDS = const_cast< OGRDataSourceType * >( this->GetOutput() );

	//Reconstruct all polygons
	std::cout << "-------------- RECONSTRUCT ALL POLYGONS FOR " << MPIConfig::Instance()->GetMyRank() <<"----------" << std::endl;
	ReconstructAllPolygons();

	//Clean Layer
	std::cout << "------------- CLEAN UNVALID FEATURES FOR " << MPIConfig::Instance()->GetMyRank() <<"----------" << std::endl;
	CleanLayer();

	//Remove border polygons from both layer
	std::cout << "------------- REMOVE OUTSIDERS FOR " << MPIConfig::Instance()->GetMyRank() <<"----------" << std::endl;
	OGRLayerType reconstructedLayer = outputOGRDS->GetLayer(reconstructedLayerName);
	OGRLayerType unvalidLayer	    = outputOGRDS->GetLayer(unvalidLayerName);
	RemoveTileBorderPolygons(reconstructedLayer);
	RemoveTileBorderPolygons(unvalidLayer);

}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::InitializeOutputDS()
{
	//We create 3 layers: one for valid polygon, one for unvalid polygon in order to repear these polygons
	OGRDataSourceType::Pointer ogrDS = OGRDataSourceType::New();
	OGRLayerType polyLayer 			  = ogrDS->CreateLayer(reconstructedLayerName, nullptr, wkbMultiPolygon);
	OGRLayerType unvalidPolygonsLayer = ogrDS->CreateLayer(unvalidLayerName, ITK_NULLPTR, wkbUnknown);

	//Create fields
	VectorOperations::CreateNewField(polyLayer, startingCoordsFieldName, OFTInteger64);
	VectorOperations::CreateNewField(polyLayer, adjStartingCoordsFieldName, OFTInteger64List);

	VectorOperations::CreateNewField(unvalidPolygonsLayer, startingCoordsFieldName, OFTInteger64);
	VectorOperations::CreateNewField(unvalidPolygonsLayer, adjStartingCoordsFieldName, OFTInteger64List);
	//Set output
	this->SetNthOutput(0, ogrDS);
}
template <class TSimplifyFunc>
typename SimplifyVectorFilter<TSimplifyFunc>::OGRFeatureType*
SimplifyVectorFilter<TSimplifyFunc>
::CreateBoundingBoxFeature()
{
	std::cout << "Create Bounding Box Polygon" << std::endl;
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );
	double ulx, uly, lrx, lry;
	//Maybe set false for bounding box
	ogrDS->GetGlobalExtent(ulx, uly, lrx, lry, true) ;

//	//Get tile information
//	int originX = std::stoi(ogrDS->ogr().GetMetadataItem("OriginTileX"));
//	int originY = std::stoi(ogrDS->ogr().GetMetadataItem("OriginTileY"));
//	int sizeX = std::stoi(ogrDS->ogr().GetMetadataItem("TileSizeX"));
//	int sizeY = std::stoi(ogrDS->ogr().GetMetadataItem("TileSizeY"));
//
//	//Modify tile origin in order to take the image origin
//	ulx = originX + ulx;
//	uly = originY + uly;
//	lrx = ulx + sizeX;
//	lry = uly + sizeY;

	//Create a feature
	OGRFeatureDefn* bbDef = new OGRFeatureDefn();
	bbDef->AddFieldDefn(new OGRFieldDefn(startingCoordsFieldName.c_str(), OFTInteger64));
	OGRFeatureType* bbFeature = new OGRFeatureType(*bbDef);
	bbFeature->ogr().SetField(startingCoordsFieldName.c_str(), -1);

	//Create polygon
	OGRLineString* bbGeom = new OGRLineString();
	bbGeom->addPoint(ulx, uly);
	bbGeom->addPoint(ulx, lry);
	bbGeom->addPoint(lrx, lry);
	bbGeom->addPoint(lrx, uly);
	bbGeom->addPoint(ulx, uly);
	//Add geomtry to feature
	bbFeature->SetGeometry(bbGeom);
	std::cout << bbGeom->exportToGML() << std::endl;
	//Remove bbPoly because internal copy is made
	delete bbGeom;

	return bbFeature;
}

template <class TSimplifyFunc>
std::vector<OGRGeometry*>
SimplifyVectorFilter<TSimplifyFunc>
::IntersectWithBoundaries(OGRFeatureType feature)
{
	//TODO:  Maybe factorize this method
	std::vector<OGRGeometry*> allEdges;

	OGRGeometry* geomRef = feature.ogr().GetGeometryRef();
	OGRGeometry* geomBB  = m_BBFeature->ogr().GetGeometryRef();

	OGRGeometry* intersectedLine = nullptr;

	std::vector<OGRGeometry*> bbEdges = VectorOperations::IntersectGeoms(geomRef, geomBB);

	//Add geometries to vector
	allEdges.insert(allEdges.end(), bbEdges.begin(), bbEdges.end());

	//Intersect with nodata features
	for(unsigned int fId = 0; fId < m_NodataLayer.GetFeatureCount(true); ++fId)
	{
		OGRFeatureType curFeature = m_NodataLayer.GetFeature(fId);
		OGRGeometry* nodataGeom = curFeature.ogr().GetGeometryRef();

		std::vector<OGRGeometry*> edges =  VectorOperations::IntersectGeoms(geomRef, nodataGeom);
		allEdges.insert(allEdges.end(), edges.begin(), edges.end());
	}

	return allEdges;
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::ConvertToEdges(OGRFeatureType feature, std::vector<OGRFeatureType> adjFeatures, OGRLayerType& insidePolygons)
{
	/**TODO: S'assurer que l'on obtient que des lines strings
	 * Si ce n'est pas le cas, on sépare de nouveau*/

	//std::cout << "Convert to edges for " << m_CurrentCoords << std::endl;
	//std::cout << "GEOM REF NOT CAST> " << feature.ogr().GetGeometryRef()->exportToGML() << std::endl;
	bool refIsParallelogram = false;
	bool adjIsParallelogram = false;
	bool isSimplify = true;
	//Geometry ref
	OGRPolygon* geomRef = VectorOperations::CastToPolygon(feature.ogr().GetGeometryRef());

	//if not a polygon, just return empty edges
	if(geomRef == nullptr)
	{
		return;
	}
	//Check if ref is rectangle
	refIsParallelogram = VectorOperations::IsParallelogram(geomRef);

	//Loop accross adj features
	for(unsigned int fId = 0; fId < adjFeatures.size(); ++fId)
	{

		//std::cout << "ID " << fId + 1 << "/" << adjFeatures.size() << std::endl;
		//Adj feature
		OGRFeatureType adjFeature = adjFeatures[fId];

		//Coordinate adj feature
		double adjCoords = adjFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());

		//Adj geom
		OGRPolygon* geomAdj = VectorOperations::CastToPolygon(adjFeature.ogr().GetGeometryRef());

		//If feature is not valid, it means we already used it
		if(VectorOperations::isFeatureValid(adjFeature) && geomAdj != nullptr)
		{
			//Check if adjacent is rectangle
			adjIsParallelogram = VectorOperations::IsParallelogram(geomAdj);

			//If one of geometries is a rectangle, we do not simplify edge
			if(refIsParallelogram || adjIsParallelogram)
			{
				isSimplify = false;
			}
			else
			{
				isSimplify = true;
			}
			//Reconstruct 2 polygons with only exterior ring to test if one is inside another one
			OGRLinearRing* adjExtRing = geomAdj->getExteriorRing();
			OGRLinearRing* refExtRing = geomRef->getExteriorRing();

			OGRPolygon* refTmp = new OGRPolygon();
			refTmp->addRing(refExtRing);

			OGRPolygon* adjTmp = new OGRPolygon();
			adjTmp->addRing(adjExtRing);

			if(adjTmp->Contains(refTmp))
			{
				//Indicate this polygon is englobed
				m_IsEnglobedMap[m_CurrentCoords] = true;
				CreateInteriorEdge(adjFeature, feature, isSimplify);
			}
			else if(refTmp->Contains(adjTmp))
			{
				m_IsEnglobedMap[adjCoords] = true;
				CreateInteriorEdge(feature, adjFeature, isSimplify);
			}
			else
			{
				//if one of feature is englobed, its adjacent is englobed too
				//(we already check if the adjacent was the englobing one)

				//Check if both are englobed, if it is, then we do not compute intersections
				if(m_IsEnglobedMap[adjCoords] ||
				   m_IsEnglobedMap[m_CurrentCoords])
				{
					//Do nothing
					std::cout << "Coord " << m_CurrentCoords << " and " << adjCoords << " both englobed" << std::endl;
				}
				else
				{
					//Compute intersected features
					IntersectFeatures(feature, adjFeature, isSimplify);
				}

			}//end if within

			//Free polygon
			delete refTmp, adjTmp;
		}//end is valid
	}//end for

}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::ConvertGeometryCollectionToEdge(OGRGeometry* intersectedLine, OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify)
{
	//Loop all geometries.
	OGRGeometryCollection* geomCollec = (OGRGeometryCollection*) intersectedLine;
	unsigned int nbGeoms = geomCollec->getNumGeometries();
	for(unsigned int geomId = 0; geomId < nbGeoms; ++geomId)
	{
		//Add edge as feature
		AddEdge(geomCollec->getGeometryRef(geomId), refFeature, adjFeature, isSimplify);
	}
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::AddEdge(OGRGeometry* intersectedLine, OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify)
{

	//Simplify
	OGRGeometry* simplifiedGeom = nullptr;
	OGRFeatureType* simplifiedFeature = nullptr;

	if(isSimplify)
	{
		m_simplifyFunc->SetInputGeom(intersectedLine);
		m_simplifyFunc->SimplifyLine();
		simplifiedGeom = m_simplifyFunc->GetOutputGeom();
	}
	else
	{
		simplifiedGeom = intersectedLine->clone();
	}

	if(simplifiedGeom != nullptr)
	{
		//Create new feature with right fields
		simplifiedFeature = new OGRFeatureType(*m_EdgeFeatureDefn);
		simplifiedFeature->SetGeometry(simplifiedGeom);

		//Set field
		double* polygonsId = new double[2];
		polygonsId[0] = refFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
		polygonsId[1] = adjFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());

		simplifiedFeature->ogr().SetField(polygonEdgeFieldName.c_str(), 2, polygonsId);

		//VectorOperations::DisplayFeature(*simplifiedFeature, false);
		//Create feature
		m_EdgeLayer.CreateFeature(*simplifiedFeature);

		//Get the updated feature
		OGRFeatureType lastFeature = (m_EdgeLayer.GetFeature(m_EdgeLayer.GetFeatureCount(true) - 1));

		//Add this feature to the map for reconstruction
		m_PolygonEdges[polygonsId[0]].push_back(lastFeature);
		m_PolygonEdges[polygonsId[1]].push_back(lastFeature);

		//Maybe free memory?
		free(polygonsId);

		//Remove simplified
		delete simplifiedGeom;

		//Delete simplifiedFeature
		delete simplifiedFeature;

	}//End if simplifiedGeom !nullptr
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::CreateInteriorEdge(OGRFeatureType englobingFeature, OGRFeatureType englobedFeature, bool isSimplify)
{
	//Add this edge to interior edges
	double englobingId = englobingFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
	double englobedId  = englobedFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());

	//Exterior ring
	OGRPolygon* englobedPolygon = (OGRPolygon*) (englobedFeature.ogr().GetGeometryRef());
	OGRLinearRing* extRing = englobedPolygon->getExteriorRing();

	//Simplify this edge
	OGRGeometry* simplified = nullptr;
	if(isSimplify)
	{
		m_simplifyFunc->SetInputGeom(extRing);
		m_simplifyFunc->SimplifyLine();
		simplified = m_simplifyFunc->GetOutputGeom();
	}
	else
	{
		simplified = extRing->clone();
	}

	//Add this interior edge to the map
	m_InteriorEdges[englobingId].push_back(simplified);

	//Keep track of who is inside who
	std::cout << "Adding " << englobedId << " in " << englobingId << std::endl;
	m_EnglobedId[englobingId].push_back(englobedId);

	//Create new feature for the interior polygon
	OGRFeatureDefn* edgeFeatureDef = new OGRFeatureDefn();
	OGRFieldDefn* fieldDef = new OGRFieldDefn(polygonEdgeFieldName.c_str(), OFTInteger64List);
	edgeFeatureDef->AddFieldDefn(fieldDef);
	OGRFeatureType simplifiedFeature(*edgeFeatureDef);
	simplifiedFeature.SetGeometry(simplified);

	//Set field
	double* polygonsId = new double[2];
	polygonsId[0] = englobingId;
	polygonsId[1] = englobedId;
	simplifiedFeature.ogr().SetField(polygonEdgeFieldName.c_str(), 2, polygonsId);
	m_PolygonEdges[englobedId].push_back(simplifiedFeature);

	delete polygonsId;
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::IntersectFeatures(OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify)
{
	//Geometry ref
	OGRGeometry* geomRef = refFeature.ogr().GetGeometryRef();

	//Adj geom
	OGRGeometry* geomAdj = adjFeature.ogr().GetGeometryRef();

	//If not inside, we intersect both geometries, and create a new one
	OGRGeometry* intersectedLine = nullptr;

	if(geomRef->Touches(geomAdj))
	{
		intersectedLine = geomRef->Intersection(geomAdj);
		OGRwkbGeometryType geomType = intersectedLine->getGeometryType();
		switch(geomType)
		{
		case wkbGeometryCollection:
			//std::cout << "geom collec" << std::endl;
			ConvertGeometryCollectionToEdge(intersectedLine,refFeature, adjFeature, isSimplify);
			break;

		default:
			//Check if we close the geometry it becomes a rectangle, then we do not simplify
			AddEdge(intersectedLine, refFeature, adjFeature, isSimplify);
			//std::cout << "add edge" << std::endl;
			break;
		}

		//Free intersectedLine
		delete intersectedLine;

	}//End if touch
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::ReconstructAllPolygons()
{
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetOutput() );
	OGRLayerType polyLayer 			  = ogrDS->GetLayer(reconstructedLayerName);
	OGRLayerType unvalidPolygonsLayer = ogrDS->GetLayer(unvalidLayerName);

	std::cout << "Number of fields in polylayer = " << polyLayer.GetLayerDefn().GetFieldCount() << std::endl;
	std::cout << "Reference count = " << polyLayer.GetLayerDefn().GetReferenceCount() << std::endl;
	std::cout << "Name = " << polyLayer.GetLayerDefn().GetFieldDefn(0)->GetNameRef() << std::endl;
	//Loop accross all features
	for (auto& edgeIt : m_PolygonEdges)
	{
		std::cout << "Reconstruction of polygon " << edgeIt.first << std::endl;
		//Construct polygon
		std::vector<double> adjCoords;
		OGRPolygon* reconstructedPolygon = ReconstructPolygon(edgeIt.first, adjCoords);

		if(reconstructedPolygon->IsValid())
		{
			//Add feature
			AddPolygon(reconstructedPolygon, edgeIt.first, adjCoords, polyLayer);
		}
		else
		{
			AddPolygon(reconstructedPolygon, edgeIt.first, adjCoords, unvalidPolygonsLayer);
		}

		//Remove polygon
		delete reconstructedPolygon;
	}
}


template <class TSimplifyFunc>
OGRPolygon*
SimplifyVectorFilter<TSimplifyFunc>
::ReconstructPolygon(double startCoords, std::vector<double>& adjCoords)
 {
	//Create vector of geometries
	//And a vector of adjCoords to update
	std::vector<OGRLineString*> unsortedGeoms = ConvertToGeometries(startCoords, adjCoords);

	//Sort vector of line string in order to create a close polygon
	std::vector<OGRLineString*> sortedGeoms = VectorOperations::SortLinesString(unsortedGeoms);

	OGRPolygon* geomPoly     = new OGRPolygon();
	OGRLinearRing* extRing   = new OGRLinearRing();

	//Loop all features
	for(unsigned int fId = 0; fId < sortedGeoms.size(); ++fId)
	{
		OGRLineString* curLine = sortedGeoms[fId];
		for(unsigned int pId = 0; pId < curLine->getNumPoints() ; ++pId)
		{
			OGRPoint curPoint;
			curLine->getPoint(pId, &curPoint);
			//extRing->addPoint(curPoint.getX(), curPoint.getY());
			extRing->addPoint(&curPoint);
		}
	}

	//Close the ring
	extRing->closeRings();

	//Add to polygon
	geomPoly->addRingDirectly(extRing);

	//Add interior rings
	std::vector<OGRLinearRing*> interiorRings = CreateInteriorRings(startCoords);
	for(unsigned int ringId = 0; ringId < interiorRings.size(); ++ringId)
	{
		geomPoly->addRingDirectly(interiorRings[ringId]);
	}


	return geomPoly;
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::AddPolygon(OGRPolygon* reconstructedPolygon, double startCoords,
			    std::vector<double> adjCoords, OGRLayerType& ogrLayer)
{
	//Add englobed id
	 std::vector<double> adjCoordsPoly = adjCoords;
	 adjCoordsPoly.insert(adjCoordsPoly.end(), m_EnglobedId[startCoords].begin(), m_EnglobedId[startCoords].end());

	//Copy adj coord to array
	double arr[adjCoordsPoly.size()];
	std::copy(adjCoordsPoly.begin(), adjCoordsPoly.end(), arr);

	//Create a feature
	OGRFeatureType polygonFeature(ogrLayer.GetLayerDefn());
	ogrLayer.CreateFeature(polygonFeature);

	polygonFeature.ogr().SetField(startingCoordsFieldName.c_str()   , startCoords);
	polygonFeature.ogr().SetField(adjStartingCoordsFieldName.c_str(), adjCoords.size(), arr);

	//Set geometry
	polygonFeature.SetGeometry(reconstructedPolygon);

	//Display
	//VectorOperations::DisplayFeature(polygonFeature, false);

	//Display Feature

	//Update
	ogrLayer.SetFeature(polygonFeature);

	//Set fields
//	std::cout << "Before = " << ogrLayer.GetLayerDefn().GetReferenceCount() << std::endl;
//	OGRFeatureType* polygonFeature = new OGRFeatureType(ogrLayer.GetLayerDefn());
//	std::cout << "After = " << ogrLayer.GetLayerDefn().GetReferenceCount() << std::endl;
//	polygonFeature->ogr().SetField(startingCoordsFieldName.c_str()   , startCoords);
//	polygonFeature->ogr().SetField(adjStartingCoordsFieldName.c_str(), adjCoords.size(), arr);
//
//
//	//Add geom to feature
//	polygonFeature->SetGeometry(reconstructedPolygon);


}

template <class TSimplifyFunc>
std::vector<OGRLinearRing*>
SimplifyVectorFilter<TSimplifyFunc>
::CreateInteriorRings(double startCoords)
{
	std::vector<OGRLinearRing*> interiorRings;

	std::vector<OGRGeometry*> interiorEdges = m_InteriorEdges[startCoords];
	for(unsigned int i = 0; i < interiorEdges.size(); ++i)
	{
		OGRLinearRing* intRing = new OGRLinearRing();
		//We clone this geometry to keep original edges
		OGRGeometry* curGeom = interiorEdges[i]->clone();
		if(curGeom->getGeometryType() == wkbLineString)
		{
			OGRLineString* lineString = (OGRLineString*) curGeom;
			for(unsigned int pId = 0; pId < lineString->getNumPoints() ; ++pId)
			{
				OGRPoint curPoint;
				lineString->getPoint(pId, &curPoint);
				//extRing->addPoint(curPoint.getX(), curPoint.getY());
				intRing->addPoint(&curPoint);
			}

			//Close ring
			intRing->closeRings();

			//Add to vector
			interiorRings.push_back(intRing);
		}
	}

	return interiorRings;
}
template <class TSimplifyFunc>
std::vector<OGRLineString*>
SimplifyVectorFilter<TSimplifyFunc>
::ConvertToGeometries(double startCoords, std::vector<double>& adjCoords)
{

	std::vector<OGRFeatureType> unsortedFeatures = m_PolygonEdges[startCoords];

	//Create vector of geometries
	std::vector<OGRLineString*> unsortedGeoms;

	//Convert to linestring
	for(unsigned int geomId = 0; geomId < m_bbEdges[startCoords].size(); ++geomId)
	{
		if(m_bbEdges[startCoords][geomId]->getGeometryType() == wkbLineString)
		{
			unsortedGeoms.push_back((OGRLineString*)m_bbEdges[startCoords][geomId]);
		}
	}

	for(unsigned int k = 0; k < unsortedFeatures.size(); ++k)
	{
		OGRFeatureType curFeature = unsortedFeatures[k];
		OGRGeometry* curGeom = curFeature.ogr().GetGeometryRef();
		OGRwkbGeometryType geomType = curGeom->getGeometryType();

		//Get coord of adjacent feature
		int nbValues = 0;
		const GIntBig* coords = curFeature.ogr().GetFieldAsInteger64List(polygonEdgeFieldName.c_str(), &nbValues);
		for(int i = 0; i < nbValues; i++)
		{
			if(coords[i] != startCoords)
			{
				//std::cout << "ADJ COORDS = " << coords[i] << std::endl;
				adjCoords.push_back(coords[i]);
			}
		}

		switch(geomType)
		{
			case wkbMultiLineString:
			{
				OGRMultiLineString* multilines = (OGRMultiLineString*)curGeom;
				for(unsigned int geomId = 0; geomId < multilines->getNumGeometries(); ++geomId)
				{
					//Theorically, we should have only linestrings
					if(multilines->getGeometryRef(geomId)->getGeometryType() == wkbLineString)
					{
						unsortedGeoms.push_back((OGRLineString*)multilines->getGeometryRef(geomId));
					}
					else
					{
						std::cout << "Not a linestring..." << std::endl;
					}
				}
				break;
			}
			case wkbLineString:
			{
				if(unsortedFeatures[k].ogr().GetGeometryRef()->getGeometryType() == wkbLineString)
				{
					unsortedGeoms.push_back((OGRLineString*)unsortedFeatures[k].ogr().GetGeometryRef());
				}
				else
				{
					std::cout << "Not a linestring..." << std::endl;
				}

				break;
			}
			default://Add nothing if it is not a linestring
				break;
		}
	}

	return unsortedGeoms;
}

template <class TSimplifyFunc>
OGRGeometry* SimplifyVectorFilter<TSimplifyFunc>
::AddInteriorRings(OGRGeometry* fixedPolygon, double startCoords)
{
	std::cout << "============= ADD INTERIOR RINGS ====================" << std::endl;
	//Get interior rings
	std::vector<OGRLinearRing*> interiorRings = CreateInteriorRings(startCoords);

	if(interiorRings.size() == 0)
	{
		std::cout << "Return copy of fixed polygon" << std::endl;
		return fixedPolygon->clone();
	}

	OGRGeometryCollection* geomCollection = new OGRGeometryCollection();

	for(unsigned int ringId = 0; ringId < interiorRings.size(); ++ringId)
	{
		std::cout << "Adding ring " << ringId << std::endl;
		//Create a polygon corresponding to this interior ring
		OGRPolygon* intPoly = new OGRPolygon();
		intPoly->addRingDirectly(interiorRings[ringId]);

		//Add polygon to geometry collection
		geomCollection->addGeometryDirectly(intPoly);

	}

	std::cout << "Display geom collection" << std::endl;
	std::cout << "Geom collection = " << geomCollection->exportToGML() << std::endl;

	//Compute difference
	OGRGeometry* updatedPolygon = fixedPolygon->Difference(geomCollection);

	if(updatedPolygon == nullptr)
	{
		std::cout << "DIFFERENCE FAILED" << std::endl;
	}

	//Delete geomCollection
	delete geomCollection;

	if(updatedPolygon->IsValid())
	{
		std::cout << "UPDATE VALID" << std::endl;
		return updatedPolygon;
	}
	else
	{
		std::cout << "UPDATE NOT VALID" << std::endl;
		return fixedPolygon->clone();
	}


}

template <class TSimplifyFunc>
void SimplifyVectorFilter<TSimplifyFunc>
::CleanLayer()
{
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetOutput() );
	OGRLayerType unvalidLayer = ogrDS->GetLayer(unvalidLayerName);
	OGRLayerType reconstructedLayer = ogrDS->GetLayer(reconstructedLayerName);

	OGRLayer* unvalidOgrLayer = &(unvalidLayer.ogr());
	OGRLayer* reconstructedOgrLayer = &(reconstructedLayer.ogr());

	//Keep track of features to delete
	std::vector<unsigned int> featuresToDelete;

	for(unsigned int i = 0; i < unvalidOgrLayer->GetFeatureCount(true); ++i)
	{
		std::cout << "Clean feature " << i + 1 << "/" << unvalidOgrLayer->GetFeatureCount(true) << std::endl;

		OGRFeature* curFeature = unvalidOgrLayer->GetFeature(i);
		int startCoords = curFeature->GetFieldAsInteger(startingCoordsFieldName.c_str());

		//Get adjacent coords
		int nbAdjacents = 0;
		const long long int* adjCoordinates = curFeature->GetFieldAsInteger64List(adjStartingCoordsFieldName.c_str(),
																				  &nbAdjacents);
		OGRFeature* cleanFeature = nullptr;

		//Dilatation
		OGRGeometry* curGeom = curFeature->GetGeometryRef();

		if(curGeom->getGeometryType() == wkbPolygon)
		{
			OGRPolygon* polygon = (OGRPolygon*) curGeom;

			//Create a self intersecting geometry
			OGRGeometry* extendedGeom = VectorOperations::ConvertToMultiPolygons(polygon);

			//Create feature
			cleanFeature = OGRFeature::CreateFeature(unvalidOgrLayer->GetLayerDefn());
			cleanFeature->SetField(startingCoordsFieldName.c_str(), startCoords);
			cleanFeature->SetField(adjStartingCoordsFieldName.c_str(), nbAdjacents, adjCoordinates);

			if(!extendedGeom->IsValid())
			{
				std::cout << "Extended not valid = " << extendedGeom->exportToGML() << std::endl;
				OGRGeometry * bufferGeom = extendedGeom->Buffer(0.0);

				if(bufferGeom == nullptr)
				{
					std::cerr << "ERROR while buffering " << polygon->exportToGML() << std::endl;
					break;
				}

				if(bufferGeom->IsValid())
				{
					//Adding interior rings if exists
					OGRGeometry* updatedGeom = AddInteriorRings(bufferGeom, startCoords);

					//Set new geometry
					cleanFeature->SetGeometryDirectly(updatedGeom);

					//Keep track of this feature to delete it after
					featuresToDelete.push_back(i);

					//Create feature for reconstructed
					reconstructedOgrLayer->CreateFeature(cleanFeature);
				}

				//Free memory
				delete bufferGeom;
			}
			else
			{
				//Adding interior rings if exists
				OGRGeometry* updatedGeom = AddInteriorRings(extendedGeom, startCoords);

				//Set new geometry
				cleanFeature->SetGeometryDirectly(updatedGeom);

				//Keep track of this feature to delete it after
				featuresToDelete.push_back(i);

				std::cout << "Before : " << updatedGeom->exportToGML() <<std::endl;
				//Create feature for reconstructed
				reconstructedOgrLayer->CreateFeature(cleanFeature);

				std::cout << "Get last feature "
						  << reconstructedOgrLayer->GetFeature(reconstructedOgrLayer->GetFeatureCount(true) - 1)
						      ->GetGeometryRef()->exportToGML() << std::endl;
			}

			//Free memory
			delete extendedGeom;

			//Delete duplicated feature
			OGRFeature::DestroyFeature(cleanFeature);

			//Delete adj coord
			delete adjCoordinates;
		}
	}//end loop features


	//Number of reconstructed features
	std::cout << "Reconstructed features for " << MPIConfig::Instance()->GetMyRank() << " = "
			  << reconstructedOgrLayer->GetFeatureCount(true) << std::endl;
	//Remove overlapping polygons
	//Indeed, when we reconstruct polygons, some can overlapp with valid polygons
	RemoveOverlappingFeatures(reconstructedLayer, featuresToDelete.size());

	std::cout << "REMOVING FIXED FEATURES FOR "<< MPIConfig::Instance()->GetMyRank() << std::endl;
	//Remove updated valid features from the unvalid layer
	for(unsigned int k = 0; k < featuresToDelete.size(); k++)
	{
		std::cout << "DELETING FEATURE : " << k +1 << "/" << featuresToDelete.size() << std::endl;
		unvalidLayer.DeleteFeature(featuresToDelete[k]);
	}

	//Check emptyness
	for(unsigned int k = 0; k  < reconstructedOgrLayer->GetFeatureCount(true); k++)
	{
		OGRGeometry* geom =  reconstructedOgrLayer->GetFeature(k)->GetGeometryRef();
		OGRPolygon* poly = (OGRPolygon*) geom;
		if(geom->IsEmpty())
		{
			std::cout << "Feature " << k << " is empty" << std::endl;
		}
	}
}

template <class TSimplifyFunc>
void SimplifyVectorFilter<TSimplifyFunc>
::RemoveOverlappingFeatures(OGRLayerType& reconstructedLayer, unsigned int nbFixedFeatures)
{
	//OGR Layer
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetOutput() );
	OGRLayer* reconstructedOgrLayer = &(reconstructedLayer.ogr());
	unsigned int nbFeatures = reconstructedOgrLayer->GetFeatureCount(true);

	// original valid feautre: 0 to (a - 1), we fixed n features: a to (a + n - 1)
	//==> a = a + n (nb features) - n (nb fixed features)
	unsigned int startId = nbFeatures - nbFixedFeatures;
	startId = 0;

	for(unsigned int fId = startId; fId < nbFeatures; ++fId)
	{
		//Current feature
		OGRFeatureType curFeature = reconstructedLayer.GetFeature(fId);
		const OGRGeometry* curGeom = curFeature.GetGeometry();
		int startCoords = curFeature.ogr().GetFieldAsInteger(startingCoordsFieldName.c_str());

		//Get adjacent feature
		std::vector<OGRFeatureType> adjFeatures = VectorOperations::GetAdjacentFeatures(ogrDS, reconstructedLayer.GetName(),
																						curFeature);

		//Loop each adj features
		for(unsigned int adjId = 0; adjId < adjFeatures.size(); adjId++)
		{
			const OGRGeometry* adjGeom = adjFeatures[adjId].GetGeometry();

			OGRGeometry* intersectedGeom = adjGeom->Intersection(curGeom);
			if(intersectedGeom->getGeometryType() == wkbPolygon)
			{
				std::cout << "Intersection is a polygon" << std::endl;
				std::cout << intersectedGeom->exportToGML() <<std::endl;
				OGRGeometry* diffGeom = curGeom->SymDifference(adjGeom);
				if(diffGeom != nullptr)
				{
					//Set new geometry for the feature
					curFeature.SetGeometry(diffGeom);

					//Update feature
					reconstructedLayer.SetFeature(curFeature);

					//Free
					delete diffGeom;
				}
			}
//			if(intersectedGeom->getGeometryType() == wkbPolygon ||
//			   intersectedGeom->getGeometryType() == wkbMultiPolygon)
//			{
			//std::cout << "DOING SYMMETRIC DIFFERENCE" << std::endl;
			//Remove adj geom from cur geom

		//	}

		}
	}
}

template <class TSimplifyFunc>
void SimplifyVectorFilter<TSimplifyFunc>
::RemoveTileBorderPolygons(OGRLayerType& layer)
{
	std::cout << "REMOVE TILE BORDER POLYGONS FOR " << layer.GetName() << std::endl;
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetOutput() );

	//Create tile geom
	OGRPolygon* tileGeom = CreateTilePolygon();

	std::vector<unsigned int> outsideFeaturesId;

	//Loop each features and check if it touches tile polygon
	for(unsigned int fid = 0; fid < layer.GetFeatureCount(true); ++fid)
	{
		OGRFeatureType currentFeature = layer.GetFeature(fid);

		const OGRGeometry* curGeom = currentFeature.GetGeometry();
		if(!curGeom->Intersects(tileGeom))
		{
			//If there is no intersection, then the feature is outside the tile
			outsideFeaturesId.push_back(fid);
		}
	}

	std::cout <<" REMOVING "<< outsideFeaturesId.size() << " OUTSIDE FEATURES" << std::endl;
	//Remove these features
	for(unsigned int fid = 0; fid < outsideFeaturesId.size(); ++fid)
	{
		layer.DeleteFeature(outsideFeaturesId[fid]);
	}

	delete tileGeom;
}

template <class TSimplifyFunc>
OGRPolygon* SimplifyVectorFilter<TSimplifyFunc>
::CreateTilePolygon()
{
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );
	int originX = std::stoi(ogrDS->ogr().GetMetadataItem("OriginTileX"));
	int originY = std::stoi(ogrDS->ogr().GetMetadataItem("OriginTileY"));
	int sizeX = std::stoi(ogrDS->ogr().GetMetadataItem("TileSizeX"));
	int sizeY = std::stoi(ogrDS->ogr().GetMetadataItem("TileSizeY"));

	//Create a feature
	OGRFeatureDefn* bbDef = new OGRFeatureDefn();
	bbDef->AddFieldDefn(new OGRFieldDefn(startingCoordsFieldName.c_str(), OFTInteger64));
	OGRFeatureType* bbFeature = new OGRFeatureType(*bbDef);
	bbFeature->ogr().SetField(startingCoordsFieldName.c_str(), -1);

	//Create polygon
	OGRLinearRing* bbGeom = new OGRLinearRing();
	bbGeom->addPoint(originX, originY);
	bbGeom->addPoint(originX, originY + sizeY);
	bbGeom->addPoint(originX + sizeX, originY + sizeY);
	bbGeom->addPoint(originX + sizeX, originY);
	bbGeom->addPoint(originX, originY);

	OGRPolygon* tilePolygon = new OGRPolygon();
	tilePolygon->addRingDirectly(bbGeom);

	return tilePolygon;
}

} //End namespace obia
} // end namespace otb

#endif
