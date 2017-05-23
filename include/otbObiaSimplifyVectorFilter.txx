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
::SimplifyVectorFilter() : m_FieldName("DN"), m_EdgeLayer(nullptr, true)
{
   this->SetNumberOfRequiredInputs(1);
   this->SetNumberOfRequiredOutputs(1);

   GDALAllRegister();

   this->ProcessObject::SetNthOutput(0, this->MakeOutput(0) );
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
	std::cout << " --------------- SIMPLIFY VECTOR FILTER ------------------" << std::endl;
	//Get the input
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );

	//Initialize output
	InitializeOutputDS();

	//Vector containing all single lines used to reconstruct polygons

	//Loop accross all geometry
	std::cout << "Generate Data for SimplifyVectorFilter" << std::endl;

	//Get number of features
//	OGRLayerType layer = ogrDS->GetLayerChecked(m_LayerName);
	OGRLayerType layer = ogrDS->GetLayer(m_LayerName);
	unsigned int nbFeatures = layer.GetFeatureCount(true);

	std::cout << "Nb features = " << nbFeatures << std::endl;

	//Creater bounding feature (bouding all the features)
	m_BBFeature = CreateBoundingBoxFeature();

	//add all features to a layer, and write it
	m_EdgeLayer = (ogrDS->CreateLayer("Edges"  , ITK_NULLPTR, wkbUnknown));
	VectorOperations::CreateNewField(m_EdgeLayer, polygonEdgeFieldName, OFTInteger64List);

	//Set definition used for next feature
	m_EdgeFeatureDefn = &(m_EdgeLayer.GetLayerDefn());

	//Maybe not usefull
	OGRLayerType insidePolygons	= ogrDS->CreateLayer("Inside"  , ITK_NULLPTR, wkbUnknown);
	VectorOperations::CreateNewField(insidePolygons, "Inside", OFTInteger64);
	VectorOperations::CreateNewField(insidePolygons, "Englobing", OFTInteger64);

	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		std::cout << "Feature id " << fId << std::endl;
		//Current feature
		OGRFeatureType curFeature = layer.GetFeature(fId);

		//Compute intersection with BB
		m_CurrentCoords = curFeature.ogr().GetFieldAsInteger64(startingCoordsFieldName.c_str());
		m_bbEdges[m_CurrentCoords] = IntersectWithBB(curFeature);

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

	//Reconstruct all polygons
	ReconstructAllPolygons();


	std::cout << "Write output" << std::endl;
	std::string  output_gml = "/space/USERS/isnard/tmp/mypolygonsSimplified.gml";
	WriteFile(output_gml, -1);
	std::cout << "End write" << std::endl;

}

template <class TSimplifyFunc>
GDALDriver*
SimplifyVectorFilter<TSimplifyFunc>
::initializeGDALDriver(std::string driverName)
{
	GDALDriver *poDriver;
	GDALAllRegister();
	poDriver = (GDALDriver*) GDALGetDriverByName(driverName.c_str());
	if( poDriver == NULL )
	{
		printf( "%s driver not available\n", driverName);
		exit( 1 );
	}

	return poDriver;

}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::InitializeOutputDS()
{
	//We create 3 layers: one for valid polygon, one for unvalid polygon in order to repear these polygons
	OGRDataSourceType::Pointer ogrDS = OGRDataSourceType::New();
	OGRLayerType polyLayer 			  = ogrDS->CreateLayer("Reconstructed", nullptr, wkbPolygon);
	OGRLayerType unvalidPolygonsLayer = ogrDS->CreateLayer("Unvalid", ITK_NULLPTR, wkbUnknown);

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

	std::cout << "ULx = " << ulx << " / ULy = " << uly << " / LRx = " << lrx << " / LRy = " << lry << std::endl;

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
::IntersectWithBB(OGRFeatureType feature)
{
	//TODO:  Maybe factorize this method
	std::vector<OGRGeometry*> bbEdges;

	OGRGeometry* geomRef = feature.ogr().GetGeometryRef();
	OGRGeometry* geomBB  = m_BBFeature->ogr().GetGeometryRef();

	OGRGeometry* intersectedLine = nullptr;

	if(geomRef->Touches(geomBB))
	{
		intersectedLine = geomRef->Intersection(geomBB);
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
					bbEdges.push_back(geomCollec->getGeometryRef(geomId));
				}
				break;
			}
			default:
				bbEdges.push_back(intersectedLine);
				break;
		}
	}//End if touch

	return bbEdges;
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::ConvertToEdges(OGRFeatureType feature, std::vector<OGRFeatureType> adjFeatures, OGRLayerType& insidePolygons)
{
	/**TODO: S'assurer que l'on obtient que des lines strings
	 * Si ce n'est pas le cas, on sépare de nouveau*/

	std::cout << "Convert to edges for " << m_CurrentCoords << std::endl;
	bool refIsRectangle = false;
	bool adjIsRectangle = false;
	bool isSimplify = true;
	//Geometry ref
	OGRPolygon* geomRef = VectorOperations::CastToPolygon(feature.ogr().GetGeometryRef());
	//if not a polygon, just return empty edges
	if(geomRef == nullptr)
	{
		return;
	}
	//Check if ref is rectangle
	refIsRectangle = VectorOperations::IsRectangle(geomRef);

	//Loop accross adj features
	for(unsigned int fId = 0; fId < adjFeatures.size(); ++fId)
	{

		std::cout << "ID " << fId + 1 << "/" << adjFeatures.size() << std::endl;
		//Adj feature
		OGRFeatureType adjFeature = adjFeatures[fId];

		//Adj geom
		OGRPolygon* geomAdj = VectorOperations::CastToPolygon(adjFeature.ogr().GetGeometryRef());

		//If feature is not valid, it means we already used it
		if(VectorOperations::isFeatureValid(adjFeature) && geomAdj != nullptr)
		{
			//Check if adjacent is rectangle
			adjIsRectangle = VectorOperations::IsRectangle(geomAdj);

			//If one of geometries is a rectangle, we do not simplify edge
			if(refIsRectangle || adjIsRectangle)
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
				CreateInteriorEdge(adjFeature, feature, isSimplify);
			}
			else if(refTmp->Contains(adjTmp))
			{
				CreateInteriorEdge(feature, adjFeature, isSimplify);
			}
			else
			{
				//Compute intersected features
				IntersectFeatures(feature, adjFeature, isSimplify);
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
	OGRGeometry* simplified = nullptr;
	OGRFeatureType* simplifiedFeature = nullptr;
	if(isSimplify)
	{
		m_simplifyFunc->SetInputGeom(intersectedLine);
		m_simplifyFunc->SimplifyLine();
		simplified = m_simplifyFunc->GetOutputGeom();
	}
	else
	{
		simplified = intersectedLine->clone();
	}

	if(simplified != nullptr)
	{
		//Create new feature with right fields
		simplifiedFeature = new OGRFeatureType(*m_EdgeFeatureDefn);
		simplifiedFeature->SetGeometry(simplified);

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
		delete simplified;

		//Delete simplifiedFeature
		delete simplifiedFeature;


	}//End if simplified !nullptr
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

			ConvertGeometryCollectionToEdge(intersectedLine,refFeature, adjFeature, isSimplify);
			break;

		default:
			AddEdge(intersectedLine, refFeature, adjFeature, isSimplify);
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
	OGRLayerType polyLayer 			  = ogrDS->GetLayer("Reconstructed");
	OGRLayerType unvalidPolygonsLayer = ogrDS->GetLayer("Unvalid");

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
		std::cout << "Adding interior ring " << ringId << std::endl;
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
	VectorOperations::DisplayFeature(polygonFeature, false);

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
		OGRGeometry* curGeom = interiorEdges[i];
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
					//Theorically, we should have only linestring
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
void SimplifyVectorFilter<TSimplifyFunc>
::WriteFile(std::string filepath, int layerId)
{
	//Get output
	OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetOutput() );
	//OGRDataSourceType::Pointer ogrDS = const_cast< OGRDataSourceType * >( this->GetInput() );
	GDALDriver*  poDriverGml = initializeGDALDriver("GML");
	GDALDataset* polyGDs = poDriverGml->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, NULL );


	if(layerId == -1)
	{
		std::cout << "Number of layer = " <<ogrDS->ogr().GetLayerCount() << std::endl;
		//Write all layers
		for(unsigned int layerId = 0; layerId < ogrDS->ogr().GetLayerCount(); ++layerId)
		{
			std::cout << "Writing layer " << ogrDS->GetLayer(layerId).ogr().GetName() << std::endl;
			std::cout << "Number feature = " << ogrDS->GetLayer(layerId).ogr().GetFeatureCount(true) << std::endl;
			OGRLayerType curLayer = ogrDS->GetLayer(layerId);
			/*OGRLayer* srcLayer 	  = &(ogrDS->GetLayer(layerId).ogr());
			OGRLayer* outputLayer = polyGDs->CreateLayer(srcLayer->GetName(), nullptr, wkbUnknown);
			for(unsigned int i = 0; i < srcLayer->GetFeatureCount(true); ++i)
			{
				outputLayer->CreateFeature(srcLayer->GetFeature(i));
			}
			 */
			polyGDs->CopyLayer(&(curLayer.ogr()), curLayer.GetName().c_str());
		}
	}
	else
	{
		otb::ogr::Layer curLayer = ogrDS->GetLayer(layerId);
		polyGDs->CopyLayer(&(curLayer.ogr()), curLayer.GetName().c_str());
		//TODO : Temporaire

//		OGRLayer* srcLayer 	  = &(layer.ogr());
//		OGRLayer* outputLayer = polyGDs->CreateLayer(layer.GetName().c_str(), nullptr, wkbUnknown);
//		for(unsigned int i = 0; i < srcLayer->GetFeatureCount(true); ++i)
//		{
//			outputLayer->CreateFeature(srcLayer->GetFeature(i));
//		}
	}


	//polyGDs->CopyLayer(&(layer.ogr()), layer.GetName().c_str(), nullptr);
	//Clean memory
	//Maybe see if we need to delete driver...
	GDALClose(polyGDs);
}
} //End namespace obia
} // end namespace otb

#endif
