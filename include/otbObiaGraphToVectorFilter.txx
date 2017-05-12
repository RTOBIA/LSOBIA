#ifndef otbObiaGraphToVectorFilter_txx
#define otbObiaGraphToVectorFilter_txx

#include "otbObiaGraphToVectorFilter.h"
#include "otbGdalDataTypeBridge.h"
#include "otbLabelImageToOGRDataSourceFilter.h"

#include "otbImage.h"
#include "otbObiaGraphToLabelImageFilter.h"

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
::GraphToVectorFilter() : m_FieldName("DN"), m_Use8Connected(false)
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
::SetInputMask(const InputGraphType *input)
{
  this->Superclass::SetNthInput(1, const_cast<InputGraphType *>(input));
}

template <class TInputGraph>
const typename GraphToVectorFilter<TInputGraph>
::InputGraphType *
GraphToVectorFilter<TInputGraph>
::GetInputMask(void)
{
  if (this->GetNumberOfInputs() < 2)
    {
    return ITK_NULLPTR;
    }

  return static_cast<const InputGraphType *>(this->Superclass::GetInput(1));
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

  typename InputGraphType::Pointer mask  =
    const_cast<InputGraphType *> (this->GetInputMask());
  if(!mask)
  {
   return;
  }
  // The input is necessarily the largest possible region.
  mask->SetRequestedRegionToLargestPossibleRegion();
}


template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::GenerateData(void)
{

	//Create the output layer for GDALPolygonize().
	ogr::DataSource::Pointer ogrDS = ogr::DataSource::New();

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
	//GDALDriver*  poDriver = initializeGDALDriver("GTiff");
	/*GDALDataset* poDataset = poDriver->Create(output_path.c_str(),nX, nY,
											   1, GDT_UInt32 , nullptr);*/

	GDALDriver*  poDriver = initializeGDALDriver("MEM");
	GDALDataset* poDataset = poDriver->Create("",nX, nY,
											   1, GDT_UInt32 , nullptr);
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
		rasterBand->RasterIO(GF_Write ,0, 0, nX, nY, fillHoleFilter->GetOutput()->GetBufferPointer(),
							  nX, nY, GDT_UInt32, 0, 0, nullptr);

	}

	//Call GDALPolygonize on this raster
	std::string layer_name = "test_layer";
	std::vector<std::string> pzOptions;
	otb::ogr::Layer layer = ogrDS->CreateLayer(layer_name, ITK_NULLPTR, wkbPolygon, pzOptions );

	OGRLayer* ogrLayer = nullptr;
	ogrLayer = &(layer.ogr());

	char* papszArgv[] {NULL};
	GDALPolygonize(rasterBand, nullptr, &(layer.ogr()), 1, papszArgv, nullptr, nullptr);
	std::cout << "Nombre polygones = " << layer.ogr().GetFeatureCount() << std::endl;

	//Clean the OGR Layer
	otb::ogr::Layer cleanedLayer = ogrDS->CreateLayer("cleaned_layer", ITK_NULLPTR, wkbPolygon, pzOptions );
	CleanOGRLayer(&(layer.ogr()), &(cleanedLayer.ogr()));

	//Remove old layer
	ogrDS->DeleteLayer(0);

	//Get output of the filter to create the ogrDs
	this->SetNthOutput(0, ogrDS);

	//Write into a file
	//Temporary
	std::string  output_gml = "/space/USERS/isnard/tmp/mypolygons.gml";
	WriteFile(output_gml);

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

	//On soustrait les subgeoms
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
		if(p.Equal(&selfIntersectingPoint)){
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
::CleanOGRLayer(OGRLayer* poLayer, OGRLayer* poLayer_cleaned)
{
	for(unsigned int i = 0; i < poLayer->GetFeatureCount(); ++i)
	{
		OGRGeometry* curGeom = poLayer->GetFeature(i)->GetGeometryRef();

		if(!curGeom->IsValid())
		{
			//std::cout << "Geometrie " << i << " non valide " << std::endl;
			if(curGeom->getGeometryType() == wkbPolygon){

				OGRPolygon* curPoly = (OGRPolygon*) curGeom;
				OGRPolygon* cleanedPoly = CleanSelfIntersectingPolygon(curPoly);
				OGRFeature* cleanFeature = new OGRFeature(poLayer->GetFeature(i)->GetDefnRef());
				cleanFeature->SetField("Label", poLayer->GetFeature(i)->GetFieldAsInteger("Label"));
				cleanFeature->SetGeometryDirectly(cleanedPoly);
				poLayer_cleaned->CreateFeature(cleanFeature);

			}
		}
		else
		{
			//Clone the feature
			poLayer_cleaned->CreateFeature(poLayer->GetFeature(i)->Clone());
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
		if(p.Equal(&listSelf[i]))
		{
			return true;
		}
	}

	return false;
}


template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::WriteFile(std::string filepath)
{
	//Get output
	auto ogrDS = this->GetOutput();
	otb::ogr::Layer layer = ogrDS->GetLayer(0);
	//TODO : Temporaire
	GDALDriver*  poDriverGml = initializeGDALDriver("GML");
	std::string  output_gml = "/space/USERS/isnard/tmp/mypolygons.gml";
	GDALDataset* polyGDs = poDriverGml->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, NULL );
	polyGDs->CopyLayer(&(layer.ogr()), "test_2", nullptr);

	//Clean memory
	//Maybe see if we need to delete driver...
	GDALClose(polyGDs);
}
} //End namespace obia
} // end namespace otb

#endif
