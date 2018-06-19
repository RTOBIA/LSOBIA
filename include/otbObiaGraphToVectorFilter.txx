/*
 * Copyright (C) 2005-2018 Centre National d'Etudes Spatiales (CNES)
 *
 * This file is part of Orfeo Toolbox
 *
 *     https://www.orfeo-toolbox.org/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef otbObiaGraphToVectorFilter_txx
#define otbObiaGraphToVectorFilter_txx

#include <sstream>
#include <string>

#include "otbObiaGraphToVectorFilter.h"
#include "otbGdalDataTypeBridge.h"
#include "otbLabelImageToOGRDataSourceFilter.h"
#include "otbObiaVectorOperations.h"

#include "otbObiaGraphToLabelImageFilter.h"
#include "otbObiaLabelImageToGraphFilter.h"
#include "otbObiaVectorOperations.h"

//RGB
#include "itkRGBPixel.h"
#include "itkLabelToRGBImageFilter.h"
#include "otbImageFileWriter.h"
//Const
#include "otbObiaConstExpr.h"

//gdal libraries
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "gdal_alg.h"
#include "gdal_utils.h"
#include "stdint.h" //needed for uintptr_t

#include "otbMPIConfig.h"

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
//  return static_cast< const OGRDataSourceType *>(
//              this->ProcessObject::GetOutput(0));
  return static_cast< const OGRDataSourceType*>(this->GetPrimaryOutput() );
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
	m_Graph = const_cast< InputGraphType * >( this->GetInput() );


	//Chain filter:
	/**First  : Convert Graph to Label
	 * Second : Convert Label to Polygons using GdalPolygonize*/

	/*****Graph to Label*******/
	//TODO Maybe add template to define input image type
	using GraphToLabelImageFilterType = otb::obia::GraphToLabelImageFilter<TInputGraph, LabelImageType>;

	unsigned int nX =  m_Graph->GetImageWidth();
	unsigned int nY =  m_Graph->GetImageHeight();
	m_Xshift = 0;//this->m_Graph->GetOriginX();
	m_Yshift = 0;//this->m_Graph->GetOriginY();

	std::cout << "Nombre de noeud : " << m_Graph->GetNumberOfNodes() << std::endl;
	auto graphToLabelFilter = GraphToLabelImageFilterType::New();
	graphToLabelFilter->SetInput(m_Graph);

	/**Fill hole in label image*/
	using FillholeFilterType = itk::GrayscaleFillholeImageFilter<LabelImageType,LabelImageType>;
	auto fillHoleFilter = FillholeFilterType::New();
	fillHoleFilter->SetInput(graphToLabelFilter->GetOutput());
	fillHoleFilter->Update();

	/***********************DEBUG : LABEL IMAGE IN RGB ************************/

//	using RGBPixelType = itk::RGBPixel<unsigned char>;
//	using RGBImageType = otb::Image<RGBPixelType, 2>;
//	using LabelToRGBFilterType = itk::LabelToRGBImageFilter<LabelImageType, RGBImageType>;
//	using RGBWriterType = otb::ImageFileWriter< RGBImageType >;
//
//	auto labelToRGBFilter = LabelToRGBFilterType::New();
//	auto rgbWriter = RGBWriterType::New();
//	std::stringstream os;
//	os.clear();
//	os << this->m_OutputDir << "RGB_IMAGE_" << MPIConfig::Instance()->GetMyRank() << ".tif";
//	rgbWriter->SetFileName(os.str());
//	labelToRGBFilter->SetInput(fillHoleFilter->GetOutput());
//	rgbWriter->SetInput(labelToRGBFilter->GetOutput());
//	rgbWriter->Update();

	/***********************END DEBUG : LABEL IMAGE IN RGB ************************/

	//Get LUT
	m_ReverseLut = graphToLabelFilter->GetReverseLut();

	std::cout << "FILL HOLE FILTER OK" << std::endl;
	//Transform this image into a GDALDataset which can be used for polygonization
	//Create a GDALDataset using MEM driver?
	/**TODO: Check the GDAL type compared to unsigned int type (16 or 32 bits?? Depends...)*/

	//Initialize driver
//	std::stringstream ss1;
//	ss1 << "Label_image_" << MPIConfig::Instance()->GetMyRank() << ".tif";
//	GDALDriver*  poDriver1 = VectorOperations::InitializeGDALDriver("GTiff");
//	GDALDataset* poDataset1 = poDriver1->Create(ss1.str().c_str(),nX, nY,
//											   1, GDT_UInt32 , nullptr);
//
//	//Create a GDALRasterBand
//	GDALRasterBand* rasterBandLabel = poDataset1->GetRasterBand(1);
//	if(rasterBandLabel != nullptr)
//	{
//		/*rasterBand->RasterIO(GF_Write ,0, 0,graph->GetImageWidth(), graph->GetImageHeight(), image_array,
//						graph->GetImageWidth(), graph->GetImageHeight(), GDT_UInt32, 0, 0, nullptr);*/
//		if(rasterBandLabel->RasterIO(GF_Write ,0, 0, nX, nY, fillHoleFilter->GetOutput()->GetBufferPointer(),
//									 nX, nY, GDT_UInt32, 0, 0, nullptr) != 0)
//		{
//			std::cout << "Error creating rasterIo" << std::endl;
//		}
//
//	}
//	GDALClose(poDataset1);
//	GeoTransform[0] /* top left x */
//	GeoTransform[1] /* w-e pixel resolution */
//  GeoTransform[2] /* 0 */
//	GeoTransform[3] /* top left y */
//	GeoTransform[4] /* 0 */
//	GeoTransform[5] /* n-s pixel resolution (negative value) */

	GDALDriver*  poDriver = VectorOperations::InitializeGDALDriver("MEM");
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
	//Call GDALPolygonize on this raster
	//Maybe pay attention to type of geometry which can be written into the layer...
	//MultiPolygons could be better, or wkbUknown to be sure?
	otb::ogr::Layer layer = ogrDS->CreateLayer(originalLayerName, ITK_NULLPTR, wkbPolygon );

	//Initialize fields
	InitializeAllFields(layer);

	char* papszArgv[] {NULL};
	GDALPolygonize(rasterBand, nullptr,  &(layer.ogr()), 0, papszArgv, nullptr, nullptr);
	std::cout << "Nombre polygones = " << layer.ogr().GetFeatureCount() << std::endl;

	//Clean the OGR Layer
	OGRLayerType cleanedLayer = ogrDS->CreateLayer(cleanedLayerName, ITK_NULLPTR, wkbUnknown);
        std::cout<<"DEBUG : adding nodatalayer to ogrDS"<<std::endl;
	OGRLayerType nodataLayer = ogrDS->CreateLayer(nodataLayerName , ITK_NULLPTR, wkbUnknown);

	//Initialize required fields
	InitializeAllFields(cleanedLayer);
	InitializeAllFields(nodataLayer);

	//Clean geometry
	CleanOGRLayer(layer, cleanedLayer, nodataLayer);

	//Remove old layer
	ogrDS->DeleteLayer(0);
	
	//Remove empty layers to avoid getting errors
	for (int i = 0; i < ogrDS->GetLayersCount(); i++)
	{
                if (ogrDS->GetLayer(i).GetFeatureCount(true) == 0)
                {
                    ogrDS->DeleteLayer(i);
                }
	}

	//Add all features
	CreateAllFeatures(cleanedLayer);

	//Set output of the filter to create the ogrDs
	//this->SetNthOutput(0, ogrDS);
	this->SetPrimaryOutput(ogrDS);

	/**DEBUG**/
	std::stringstream ss;
	ss << "mypolygons_" << MPIConfig::Instance()->GetMyRank() << ".gml";
	std::string  output_gml = ss.str();
	VectorOperations::WriteOGRDataSource(ogrDS, output_gml, -1);

	//Validate layer?
	std::cout << "Check layer = " << ogrDS->GetLayerChecked(cleanedLayerName).GetFeatureCount(true) << std::endl;
	
	//Clear memory
	GDALClose(poDataset);
}

//Create a field
template <class TInputGraph>
void
GraphToVectorFilter<TInputGraph>
::InitializeAllFields(OGRLayerType& poLayer)
{
	//ATTENTION: Label field must be the first in order to the polygonize function use index 0 (label field)
	//to write pixel value
	VectorOperations::CreateNewField(poLayer, labelFieldName			  , OFTInteger);
	VectorOperations::CreateNewField(poLayer, validPolygonFieldName	      , OFTString);
	VectorOperations::CreateNewField(poLayer, startingCoordsFieldName     , OFTInteger64);
	VectorOperations::CreateNewField(poLayer, adjStartingCoordsFieldName  , OFTInteger64List);
}

//
//template <class TInputGraph>
//std::vector<OGRPoint>
//GraphToVectorFilter<TInputGraph>
//::GetSelfIntersectingPoints(OGRPolygon* ogrPolygon)
//{
//	//std::cout << "Self Intersecting Points " << std::endl;
//	std::map<std::pair<double, double>, std::vector<OGRPoint>> mapPoints;
//	std::vector<OGRPoint> selfIntersectedPoints;
//	OGRCurve* ogrCurve = ogrPolygon->getExteriorRingCurve();
//	if(ogrCurve != nullptr){
//		OGRPointIterator* ptIt = ogrCurve->getPointIterator();
//		OGRPoint p;
//			while(ptIt->getNextPoint(&p)){
//				//std::cout << "Point " << p.exportToGML() << std::endl;
//				std::pair<double, double> curCoord;
//
//				curCoord.first  = p.getX();
//				curCoord.second = p.getY();
//				mapPoints[curCoord].push_back(p);
//			}
//
//			//Suppression du premier et dernier point
//			ptIt = ogrPolygon->getExteriorRingCurve()->getPointIterator();
//			ptIt->getNextPoint(&p);
//			std::pair<double, double> firstCoord;
//			firstCoord.first = p.getX();
//			firstCoord.second = p.getY();
//			mapPoints.erase(firstCoord);
//
//			//Parcours de la map
//		    for (auto& x: mapPoints) {
//		            if(x.second.size() > 1){
//		            	selfIntersectedPoints.insert(selfIntersectedPoints.end(), x.second.begin(), x.second.end());
//		            }
//		    }
//	}
//
//
//    return selfIntersectedPoints;
//}
//
//
//template <class TInputGraph>
//OGRPolygon* GraphToVectorFilter<TInputGraph>
//::CleanSelfIntersectingPolygon(OGRPolygon* ogrPolygon)
//{
//	//Nouveau polygon
//	OGRPolygon* cleanPoly = new OGRPolygon();
//	std::vector<OGRPolygon*> subGeoms;
//	OGRLinearRing* linearRing = new OGRLinearRing;
//
//	//Get number of selfs intersecting points
//	std::vector<OGRPoint> selfIntersectingPoints = GetSelfIntersectingPoints(ogrPolygon);
//
//	//std::cout << "Nombre points Self Intersecting " << selfIntersectingPoints.size() << std::endl;
//	//std::cout << "Point Self " << (*selfIntersectingPoints.begin()).exportToGML() << std::endl;
//
//	//Parcours de chaque points du polygone
//	OGRPointIterator* ptIt = ogrPolygon->getExteriorRing()->getPointIterator();
//	OGRPoint curPoint;
//	while(ptIt->getNextPoint(&curPoint))
//	{
//		//std::cout << curPoint.exportToGML() << std::endl;
//		//Ajout du point
//		linearRing->addPoint(&curPoint);
//		if(IsSelfIntersecting(curPoint, selfIntersectingPoints))
//		{
//			//std::cout << "############ Creation nouvelle geometrie" << std::endl;
//			OGRPolygon* subGeom = CreateSubGeomtry(ptIt, curPoint);
//			subGeoms.push_back(subGeom);
//			//std::cout << "{######################################" << std::endl;
//		}
//
//		//Si ce point appartient à self intersecting, alors :
//		//On ajoute le point au polygone nettoyé
//		//On crée un nouveau polygone correspondant à la forme self intersected
//		//On passe au point d'après
//
//	}
//
//	//Add ring directly
//	cleanPoly->addRingDirectly(linearRing);
//
//	//On soustrait les subgeoms (dans le cas où on a un polygone a l'intérieur d'un autre par exemple)
//	for(unsigned int k = 0; k < subGeoms.size(); ++k)
//	{
//		//TODO
//		delete subGeoms[k];
//	}
//	return cleanPoly;
//}
//
//
//template <class TInputGraph>
//OGRPolygon* GraphToVectorFilter<TInputGraph>
//::CreateSubGeomtry(OGRPointIterator* ptIt, OGRPoint selfIntersectingPoint)
//{
//	OGRPolygon* subGeom = new OGRPolygon();
//	OGRLinearRing linearRing;
//
//	linearRing.addPoint(&selfIntersectingPoint);
//
//	//Parcours des points jusqu'à retomber sur le self intersecting points
//	OGRPoint p;
//	bool isClosed = false;
//	while(!isClosed)
//	{
//		ptIt->getNextPoint(&p);
//		if(p.Equals(&selfIntersectingPoint)){
//			/*std::cout << "p " << p.getX() << "/" << p.getY()
//					  << " compared " << selfIntersectingPoint.getX() << "/" << selfIntersectingPoint.getY() << std::endl;*/
//			isClosed = true;
//		}else{
//			//std::cout << "p " << p.getX() << "/" << p.getY() << std::endl;
//			linearRing.addPoint(&p);
//		}
//	}
//
//	subGeom->addRing(&linearRing);
//	//std::cout << "Sub Geom " << subGeom->getExteriorRing()->getNumPoints() << std::endl;
//
//	return subGeom;
//	//Add points
//
//}


template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::CleanOGRLayer(OGRLayerType poLayer, OGRLayerType& poLayer_cleaned, OGRLayerType& polayer_nodata)
{
	std::cout << "================ CLEAN OGR LAYER ====================" << std::endl;
	OGRLayer* cleanedLayer  = &(poLayer_cleaned.ogr());
	OGRLayer* nodataLayer	= &(polayer_nodata.ogr());
	OGRLayer* originalLayer = &(poLayer.ogr());

	for(unsigned int i = 0; i < originalLayer->GetFeatureCount(); ++i)
	{
		OGRFeature* curFeature = originalLayer->GetFeature(i);
		OGRFeature* cleanFeature = nullptr;
		OGRGeometry* curGeom = curFeature->GetGeometryRef();
		int labelValue = curFeature->GetFieldAsInteger(labelFieldName.c_str());

		if(labelValue != 0)
		{
			if(!curGeom->IsValid())
			{

				//Dilatation
				OGRGeometry * bufferGeom = curGeom->Buffer(0.0);
				if(bufferGeom->IsValid())
				{
					cleanFeature = OGRFeature::CreateFeature(originalLayer->GetLayerDefn());
					cleanFeature->SetField(labelFieldName.c_str(), labelValue);
					OGRErr geomErr = cleanFeature->SetGeometryDirectly(bufferGeom);
					if(geomErr != 0)
					{
						std::cout << "Error setting geometry for " << i << std::endl;
					}
				}

			}
			else
			{
				cleanFeature = OGRFeature::CreateFeature(originalLayer->GetLayerDefn());
				cleanFeature->SetField(labelFieldName.c_str(), labelValue);
				OGRErr geomErr = cleanFeature->SetGeometry(curGeom);
				if(geomErr != 0)
				{
					std::cout << "Error setting geometry for " << i << std::endl;
				}
			}

			//Add this feature to the clean layer
			OGRErr err = cleanedLayer->CreateFeature(cleanFeature);
			if(err != 0)
			{
				std::cout << err << " : Error creating feature " << i << std::endl;
			}

		}
		else
		{
                        std::cout<<"labelValue == 0 : nodata"<<std::endl;
			//We need to add the nodata value feature in order to keep the boundaries
			//std::cout << "Nodata value GEOMETRY for " << MPIConfig::Instance()->GetMyRank()<< std::endl;
			//std::cout << curFeature->GetGeometryRef()->exportToGML() << std::endl;
			cleanFeature = OGRFeature::CreateFeature(originalLayer->GetLayerDefn());
			cleanFeature->SetField(labelFieldName.c_str(), -1);
			OGRErr geomErr = cleanFeature->SetGeometry(curFeature->GetGeometryRef());
			if(geomErr != 0)
			{
				std::cout << "Error setting geometry for " << i << std::endl;
			}
			//Add to nodatalayer
			OGRErr err = nodataLayer->CreateFeature(cleanFeature);
			if(err != 0)
			{
				std::cout << err << " : Error creating feature " << i << std::endl;
			}

		}

		//Destroy feature because it is not consummed by CreateFeature
	    OGRFeature::DestroyFeature(cleanFeature);
	}
}

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::CreateAllFeatures(OGRLayerType& poLayer)
{
	std::cout << "Create all features for " << MPIConfig::Instance()->GetMyRank() << std::endl;
	OGRLayer* ogrLayer = &(poLayer.ogr());
	//Loop all polygons
	unsigned int nbFeatures = poLayer.GetFeatureCount(true);
	for(unsigned int fId = 0; fId < nbFeatures; ++fId)
	{
		//Get geometry
		OGRFeature* curFeature = ogrLayer->GetFeature(fId);

		if(curFeature != nullptr)
		{
			//Check if geometry valid?
			if(curFeature->GetGeometryRef()->IsValid())
			{
				//If valid, update feature and add
				UpdateFeatureFields(curFeature);

				//Update the feature
				poLayer.SetFeature(curFeature);
			}
		}
	}

	std::cout << "End creating all features for " << MPIConfig::Instance()->GetMyRank() << std::endl;

}

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::UpdateFeatureFields(OGRFeature* curFeature)
 {
	//Get geometry
	int label = curFeature->GetFieldAsInteger(labelFieldName.c_str());

	//Label -1 is for nodata features, we do not need to update it
	if(label == - 1)
	{
		return;
	}

	int reverseLabel = m_ReverseLut[label];
	//Get associated node
	//The label image has been created using label = nodeId + 1
	//So the reverse operation is nodeId = label - 1
	auto node = m_Graph->GetNodeAt(reverseLabel - 1);

	if(node == nullptr)
	{
		std::cout << "NODE " << label - 1 << " is null" << std::endl;
		return;
	}

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
		//std::cout << "Get node at " << edgeIt->m_TargetId << std::endl;
		auto adjNode = m_Graph->GetNodeAt(edgeIt->m_TargetId);

		if(adjNode == nullptr)
		{
			std::cout << "ADJ NODE of " << startingCoords << " is null" << std::endl;
			return;
		}

		//Add to field
		double adjStartingCoords = adjNode->GetFirstPixelCoords();
		adjacentCoords[cpt] =  adjStartingCoords;
		cpt++;
	}

	//Add starting coord of each adjacent node
	curFeature->SetField(adjStartingCoordsFieldName.c_str(),  nbAdjacents, adjacentCoords);
	curFeature->SetField(validPolygonFieldName.c_str(), "true");

	//delete adjacent coords?
	delete adjacentCoords;
	adjacentCoords = nullptr;
 }

template <class TInputGraph>
void GraphToVectorFilter<TInputGraph>
::UpdateOutputInformation()
{
	Superclass::UpdateOutputInformation();

	std::cout << "Update output" << std::endl;
    // Mise a jour du graphe en fonction de son header
	auto ogrDS = const_cast<OGRDataSourceType*>(this->GetOutput());

    itk::ModifiedTimeType t1;
    t1 = this->GetMTime();
    ogrDS->SetPipelineMTime(t1);

}

} //End namespace obia
} // end namespace otb

#endif
