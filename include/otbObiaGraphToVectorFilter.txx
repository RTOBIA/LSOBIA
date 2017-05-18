#ifndef otbObiaGraphToVectorFilter_txx
#define otbObiaGraphToVectorFilter_txx

#include "otbObiaGraphToVectorFilter.h"
#include "otbGdalDataTypeBridge.h"
#include "otbLabelImageToOGRDataSourceFilter.h"

#include "otbImage.h"
#include "otbObiaGraphToLabelImageFilter.h"

//Const
#include "otbObiaConstExpr.h"

//gdal libraries
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "gdal_alg.h"
#include "gdal_utils.h"
#include "stdint.h" //needed for uintptr_t

namespace otb
{
namespace obia
{

template <class TInputGraph>
GraphToVectorFilter<TInputGraph>
::GraphToVectorFilter() : m_FieldName("DN"), m_Use8Connected(false), m_Xshift(0), m_Yshift(0)
{
   this->SetNumberOfRequiredInputs(1);
   this->SetNumberOfRequiredOutputs(1);

   GDALAllRegister();

   this->ProcessObject::SetNthOutput(0, this->MakeOutput(0) );
}


template <class TInputGraph>
typename GraphToVectorFilter<TInputGraph>::DataObjectPointer
GraphToVectorFilter<TInputGraph>
::MakeOutput(DataObjectPointerArraySizeType itkNotUsed(idx))
{
  return static_cast< DataObjectPointer >(OGRDataSourceType::New().GetPointer());
}

template <class TInputGraph>
const typename GraphToVectorFilter<TInputGraph>::OGRDataSourceType *
GraphToVectorFilter<TInputGraph>
::GetOutput()
{
  return static_cast< const OGRDataSourceType *>(
              this->ProcessObject::GetOutput(0));
}

template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::SetInput(const InputGraphType *input)
{
  this->Superclass::SetNthInput(0, const_cast<InputGraphType *>(input));
}

template <class TInputGraph>
const typename GraphToVectorFilter<TInputGraph>
::InputGraphType *
GraphToVectorFilter<TInputGraph>
::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
    {
    return ITK_NULLPTR;
    }

  return static_cast<const InputGraphType *>(this->Superclass::GetInput(0));
}

template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::GenerateInputRequestedRegion(void)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename InputGraphType::Pointer input  =
    const_cast<InputGraphType *> (this->GetInput());

  if ( !input )
    {
    return;
    }
  // The input is necessarily the largest possible region.
  input->SetRequestedRegionToLargestPossibleRegion();

}


template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::GenerateData(void)
{
	//Output
	//Use Self::Pointer in
	OGRDataSourceType::Pointer ogrDS = ogr::DataSource::New();


	//Initialize graph pointer
	m_graph = const_cast< InputGraphType * >( this->GetInput() );

	//Chain filter:
	/**First  : Convert Graph to Label
	 * Second : Convert Label to Polygons using GdalPolygonize*/

	/*****Graph to Label*******/
	//TODO Maybe add template to define input image type
    using LabelPixelType = unsigned int;
    using LabelImageType = otb::Image< LabelPixelType, 2 >;
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<TInputGraph, LabelImageType>;

	auto graph = this->GetInput();
	unsigned int nX =  graph->GetImageWidth();
	unsigned int nY =  graph->GetImageHeight();
	std::cout << "Nombre de noeud : " << graph->GetNumberOfNodes() << std::endl;
	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	graphToLabelFilter->SetInput(graph);

	/**Fill hole in label image*/
	using FillholeFilterType = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;
	auto fillHoleFilter = FillholeFilterType::New();
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	fillHoleFilter->Update();

	//Transform this image into a GDALDataset which can be used for polygonization
	//Create a GDALDataset using MEM driver?
	/**TODO: Check the GDAL type compared to unsigned int type (16 or 32 bits?? Depends...)*/

	//Initialize driver
	/*std::string  output_img = "/space/USERS/isnard/tmp/mygdal.tif";
	GDALDriver*  poDriver = initializeGDALDriver("GTiff");
	GDALDataset* poDataset = poDriver->Create(output_img.c_str(),nX, nY,
											   1, GDT_UInt32 , nullptr);*/


//	GeoTransform[0] /* top left x */
//	GeoTransform[1] /* w-e pixel resolution */
//  GeoTransform[2] /* 0 */
//	GeoTransform[3] /* top left y */
//	GeoTransform[4] /* 0 */
//	GeoTransform[5] /* n-s pixel resolution (negative value) */

	GDALDriver*  poDriver = initializeGDALDriver("MEM");
	GDALDataset* poDataset = poDriver->Create("",nX, nY,
											   1, GDT_UInt32 , nullptr);

	//Adding a geo transform to get right coordinates for polygons
	double* geotransform = new double[6];
	geotransform[0] = m_Xshift;
	geotransform[1] = 1;
	geotransform[2] = 0;
	geotransform[3] = m_Yshift;
	geotransform[4] = 0;
	geotransform[5] = 1;

	poDataset->SetGeoTransform(geotransform);

	if(poDataset == nullptr)
	{
		std::cout << "Problem with creating gdal dataset...\n" << std::endl;
		exit(1);
	}

	//Create a GDALRasterBand
	GDALRasterBand* rasterBand = poDataset->GetRasterBand(1);
	if(rasterBand != nullptr)
	{
		/*rasterBand->RasterIO(GF_Write ,0, 0,graph->GetImageWidth(), graph->GetImageHeight(), image_array,
						graph->GetImageWidth(), graph->GetImageHeight(), GDT_UInt32, 0, 0, nullptr);*/
		if(rasterBand->RasterIO(GF_Write ,0, 0, nX, nY, fillHoleFilter->GetOutput()->GetBufferPointer(),
							  nX, nY, GDT_UInt32, 0, 0, nullptr) != 0)
		{
			std::cout << "Error creating rasterIo" << std::endl;
		}

	}


	//Initialise the OGR DS

	//Create the output layer for GDALPolygonize().
	//std::cout << "Initialize OGR DS" << std::endl;
	//ogrDS = ogr::DataSource::New(poDataset, ogr::DataSource::Modes::Overwrite);
	//std::cout << "End initialise OGR DS" << std::endl;

	//Call GDALPolygonize on this raster
	//Maybe pay attention to type of geometry which can be written into the layer...
	//MultiPolygons could be better, or wkbUknown to be sure?
	std::cout << "Creating layer" << std::endl;
	otb::ogr::Layer layer = ogrDS->CreateLayer(originalLayerName, ITK_NULLPTR, wkbPolygon );

	std::cout << "End create layer" << std::endl;
	//Initialize fields
	InitializeAllFields(layer);

	char* papszArgv[] {NULL};
	GDALPolygonize(rasterBand, nullptr,  &(layer.ogr()), 0, papszArgv, nullptr, nullptr);
	std::cout << "Nombre polygones = " << layer.ogr().GetFeatureCount() << std::endl;

	//Clean the OGR Layer
	OTBLayer cleanedLayer = ogrDS->CreateLayer(cleanedLayerName, ITK_NULLPTR, wkbUnknown);

	InitializeAllFields(cleanedLayer);
	/*CreateNewField(cleanedLayer, labelFieldName				, OFTInteger);
	CreateNewField(cleanedLayer, startingCoordsFieldName	, OFTInteger64);
	CreateNewField(cleanedLayer, adjStartingCoordsFieldName	, OFTInteger64List);*/

	//Clean geometry
	CleanOGRLayer(layer, cleanedLayer);

	//Remove old layer
	ogrDS->DeleteLayer(0);

	//Add all features
	CreateAllFeatures(cleanedLayer);

	//Set output of the filter to create the ogrDs
	this->SetNthOutput(0, ogrDS);


	//Write into a file
	//Temporary
	std::string  output_gml = "/space/USERS/isnard/tmp/mypolygons.gml";
	WriteFile(output_gml);

	//Validate layer?
	std::cout << "Check layer = " << ogrDS->GetLayerChecked(cleanedLayerName).GetFeatureCount(true) << std::endl;

	//Clear memory
	GDALClose(poDataset);



}

template <class TInputGraph>
GDALDriver*
GraphToVectorFilter<TInputGraph>
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

//Create a field
template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::InitializeAllFields(OTBLayer& poLayer)
{
	//ATTENTION: Label field must be the first in order to the polygonize function use index 0 (label field)
	//to write pixel value
	CreateNewField(poLayer, labelFieldName			  , OFTInteger);
	CreateNewField(poLayer, validPolygonFieldName	  , OFTString);
	CreateNewField(poLayer, startingCoordsFieldName   , OFTInteger64);
	CreateNewField(poLayer, adjStartingCoordsFieldName, OFTInteger64List);
}

//Create a field
template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::CreateNewField(otb::ogr::Layer poLayer, std::string fieldName, OGRFieldType fieldType)
{
	OGRFieldDefn fieldDef(fieldName.c_str(), fieldType);
	poLayer.CreateField(fieldDef, true);

}
template <class TInputGraph>
std::vector<OGRPoint>
GraphToVectorFilter<TInputGraph>
::GetSelfIntersectingPoints(OGRPolygon* ogrPolygon)
{
	//std::cout << "Self Intersecting Points " << std::endl;
	std::map<std::pair<double, double>, std::vector<OGRPoint>> mapPoints;
	std::vector<OGRPoint> selfIntersectedPoints;
	OGRCurve* ogrCurve = ogrPolygon->getExteriorRingCurve();
	if(ogrCurve != nullptr){
		OGRPointIterator* ptIt = ogrCurve->getPointIterator();
		OGRPoint p;
			while(ptIt->getNextPoint(&p)){
				//std::cout << "Point " << p.exportToGML() << std::endl;
				std::pair<double, double> curCoord;

				curCoord.first  = p.getX();
				curCoord.second = p.getY();
				mapPoints[curCoord].push_back(p);
			}

			//Suppression du premier et dernier point
			ptIt = ogrPolygon->getExteriorRingCurve()->getPointIterator();
			ptIt->getNextPoint(&p);
			std::pair<double, double> firstCoord;
			firstCoord.first = p.getX();
			firstCoord.second = p.getY();
			mapPoints.erase(firstCoord);

			//Parcours de la map
		    for (auto& x: mapPoints) {
		            if(x.second.size() > 1){
		            	selfIntersectedPoints.insert(selfIntersectedPoints.end(), x.second.begin(), x.second.end());
		            }
		    }
	}


    return selfIntersectedPoints;
}


template <class TInputGraph>
OGRPolygon* GraphToVectorFilter<TInputGraph>
::CleanSelfIntersectingPolygon(OGRPolygon* ogrPolygon)
{



	//Nouveau polygon
	OGRPolygon* cleanPoly = new OGRPolygon();
	std::vector<OGRPolygon*> subGeoms;
	OGRLinearRing linearRing;

	//Get number of selfs intersecting points
	std::vector<OGRPoint> selfIntersectingPoints = GetSelfIntersectingPoints(ogrPolygon);

	//std::cout << "Nombre points Self Intersecting " << selfIntersectingPoints.size() << std::endl;
	//std::cout << "Point Self " << (*selfIntersectingPoints.begin()).exportToGML() << std::endl;

	//Parcours de chaque points du polygone
	OGRPointIterator* ptIt = ogrPolygon->getExteriorRing()->getPointIterator();
	OGRPoint curPoint;
	while(ptIt->getNextPoint(&curPoint))
	{
		//std::cout << curPoint.exportToGML() << std::endl;
		//Ajout du point
		linearRing.addPoint(&curPoint);
		if(IsSelfIntersecting(curPoint, selfIntersectingPoints))
		{
			//std::cout << "############ Creation nouvelle geometrie" << std::endl;
			OGRPolygon* subGeom = CreateSubGeomtry(ptIt, curPoint);
			subGeoms.push_back(subGeom);
			//std::cout << "{######################################" << std::endl;
		}

		//Si ce point appartient à self intersecting, alors :
		//On ajoute le point au polygone nettoyé
		//On crée un nouveau polygone correspondant à la forme self intersected
		//On passe au point d'après

	}
	cleanPoly->addRing(&linearRing);

	//On soustrait les subgeoms (dans le cas où on a un polygone a l'intérieur d'un autre par exemple)
	for(unsigned int k = 0; k < subGeoms.size(); ++k)
	{
		//TODO
	}
	return cleanPoly;
}


template <class TInputGraph>
OGRPolygon* GraphToVectorFilter<TInputGraph>
::CreateSubGeomtry(OGRPointIterator* ptIt, OGRPoint selfIntersectingPoint)
{
	OGRPolygon* subGeom = new OGRPolygon();
	OGRLinearRing linearRing;

	linearRing.addPoint(&selfIntersectingPoint);

	//Parcours des points jusqu'à retomber sur le self intersecting points
	OGRPoint p;
	bool isClosed = false;
	while(!isClosed)
	{
		ptIt->getNextPoint(&p);
		if(p.Equals(&selfIntersectingPoint)){
			/*std::cout << "p " << p.getX() << "/" << p.getY()
					  << " compared " << selfIntersectingPoint.getX() << "/" << selfIntersectingPoint.getY() << std::endl;*/
			isClosed = true;
		}else{
			//std::cout << "p " << p.getX() << "/" << p.getY() << std::endl;
			linearRing.addPoint(&p);
		}
	}

	subGeom->addRing(&linearRing);
	//std::cout << "Sub Geom " << subGeom->getExteriorRing()->getNumPoints() << std::endl;

	return subGeom;
	//Add points

}


template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::CleanOGRLayer(otb::ogr::Layer poLayer, otb::ogr::Layer& poLayer_cleaned)
{
	OGRLayer* cleanedLayer  = &(poLayer_cleaned.ogr());
	OGRLayer* originalLayer = &(poLayer.ogr());

	for(unsigned int i = 0; i < originalLayer->GetFeatureCount(); ++i)
	{
		OGRGeometry* curGeom = originalLayer->GetFeature(i)->GetGeometryRef();

		if(!curGeom->IsValid())
		{
			if(curGeom->getGeometryType() == wkbPolygon){

				OGRPolygon* curPoly = (OGRPolygon*) curGeom;
				OGRPolygon* cleanedPoly = CleanSelfIntersectingPolygon(curPoly);
				OGRFeature* cleanFeature = new OGRFeature(originalLayer->GetFeature(i)->GetDefnRef());
				cleanFeature->SetField(labelFieldName.c_str(),
										originalLayer->GetFeature(i)->GetFieldAsInteger(labelFieldName.c_str()));
				if( cleanFeature->SetGeometryDirectly(cleanedPoly) == 0)
				{
					if(cleanedLayer->CreateFeature(cleanFeature) != 0)
					{
						std::cout << "Error creating feature" << std::endl;
					}
				}
				else
				{
					std::cout << "Error setting geometry directly" << std::endl;
				}
			}
		}
		else
		{
			//Clone the feature
			if(cleanedLayer->CreateFeature(originalLayer->GetFeature(i)->Clone()) != 0)
			{
				std::cout << "Error creating feature" << std::endl;
			}
		}

	}
	std::cout << "End Clean" << std::endl;
}


template <class TInputGraph>
bool GraphToVectorFilter<TInputGraph>
::IsSelfIntersecting(OGRPoint p, std::vector<OGRPoint> listSelf)
{
	for(unsigned int i = 0; i < listSelf.size(); ++i)
	{
		if(p.Equals(&listSelf[i]))
		{
			return true;
		}
	}

	return false;
}

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::CreateAllFeatures(OTBLayer& poLayer)
{
	std::cout << "Create all features" << std::endl;

	//Loop all polygons
	unsigned int nbFeatures = poLayer.GetFeatureCount(true);
	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		//Get geometry
		OGRLayer* ogrLayer = &(poLayer.ogr());
		OGRFeature* curFeature = ogrLayer->GetFeature(fId);

		//Check if geometry valid?
		if(curFeature->GetGeometryRef()->IsValid())
		{
			//If valid, update feature and add
			UpdateFeatureFields(curFeature);

			//Update the feature
			poLayer.SetFeature(curFeature);
		}
		else
		{
			std::cout << "Feature " << fId << " not valid " << std::endl;
		}

	}

}

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::UpdateFeatureFields(OGRFeature* curFeature)
 {
	//Get geometry
	OGRGeometry* curGeom  = curFeature->GetGeometryRef();
	int label = curFeature->GetFieldAsInteger(labelFieldName.c_str());

	//std::cout << "Label " << label << std::endl;
	//std::cout << "Number of fields = " << curFeature->GetFieldCount() << std::endl;
	//Get associated node
	//The label image has been created using label = nodeId + 1
	//So the reverse operation is nodeId = label - 1
	auto node = m_graph->GetNodeAt(label - 1);

	//Add starting coord as field
	long long int startingCoords = node->GetFirstPixelCoords();
	//std::cout << "startingCoords = " << startingCoords << std::endl;
	curFeature->SetField(startingCoordsFieldName.c_str(), startingCoords);

	double* adjacentCoords = new double[node->m_Edges.size()];
	unsigned int nbAdjacents = node->m_Edges.size();
	unsigned int cpt = 0;
	//Loop accross each adjacent
	for(auto edgeIt = node->m_Edges.begin(); edgeIt != node->m_Edges.end(); edgeIt++)
	{
		//Adjacent node
		auto adjNode = m_graph->GetNodeAt(edgeIt->m_TargetId);

		//Add to field
		double adjStartingCoords = adjNode->GetFirstPixelCoords();
		adjacentCoords[cpt] =  adjStartingCoords;
		cpt++;
	}

	//Add starting coord of each adjacent node
	curFeature->SetField(adjStartingCoordsFieldName.c_str(),  nbAdjacents, adjacentCoords);
	curFeature->SetField(validPolygonFieldName.c_str(), "true");
 }

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::WriteFile(std::string filepath)
{
	//Get output
	auto ogrDS = this->GetOutput();
	otb::ogr::Layer layer = ogrDS->GetLayer(0);
	std::cout << "Layer name = " << layer.GetName() << std::endl;

	GDALDriver*  poDriverGml = initializeGDALDriver("GML");
	GDALDataset* polyGDs = poDriverGml->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, NULL );
	polyGDs->CopyLayer(&(layer.ogr()), layer.ogr().GetName(), nullptr);
	//Clean memory
	//Maybe see if we need to delete driver...
	GDALClose(polyGDs);
}
} //End namespace obia
} // end namespace otb

#endif
