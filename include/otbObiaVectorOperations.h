#ifndef otbObiaVectorOperations_h
#define otbObiaVectorOperations_h

#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <unordered_set>

#include <otbOGRDataSourceWrapper.h>

/**
\file otbObiaVectorOperations.h
\brief This file define the VectorOperations class which allows some operations on vector
 (like extract adjacent features)
*/
namespace otb
{
namespace obia
{
/** \class VectorOperations
 *	\brief Class which allow some operations of vector, and features (like converting to polygon, extracting  features...)
 *
 */
    class VectorOperations
    {

      public:

      /* Some convenient alias */
      using OGRDataSourceType		 = otb::ogr::DataSource;
      using OGRFeatureType			 = otb::ogr::Feature;
      using OGRLayerType			 = otb::ogr::Layer;

      // \brief This method extract feature from the OGRDatasource with the field value given
      //
      // @param ogrDS: OGRDataSource
      // @param layerName : layer where the feature is
      // @param field name : name of the field to look at
      // @param field value : value of the field to look at
      // @return feature: return the feature (null if not found)
      template<typename T>
      static OGRFeatureType GetFeatureAt(OGRDataSourceType* ogrDS,
    		  	  	  	  	   	   	   	  const std::string layerName,
										  const std::string fieldName,
										  const T fieldValue);

      // \brief This method extract all features which satisfy the field value given
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param field name: name of the field to look at
      // @param field value : value of the field to look at
      template<typename T>
      static std::vector<OGRFeatureType> GetAllFeaturesAt(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  const std::string layerName,
														  const std::string fieldName,
														  const T fieldValue);

      //  \brief  This method extract all adjacent feature to the given feature
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param feature : reference feature
      // @return list of adjacent features
      static std::vector<OGRFeatureType> GetAdjacentFeatures(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  	 const std::string layerName,
															 const OGRFeatureType feature);

      /** \brief Check if a feature is valid (meaning already used for interscting for example)
       * @param: input feature
       * @return : true or false*/
      static bool isFeatureValid(const OGRFeatureType feature);

      /** \brief  Sort linestring in order to reconstruct a polygon
       * @param : Unsorted linestring
       * @return : vector of sorted linestring*/
      static std::vector<OGRLineString*> SortLinesString(std::vector<OGRLineString*> unsortedGeoms);

      /** \brief Cast geometry to polygon geometry
       * @param: Geometry
       * @return: casted polygon*/
      static OGRPolygon* CastToPolygon(OGRGeometry* geomRef);

      /** \brief Create a new field to the input layer*/
      static void CreateNewField(OGRLayerType& poLayer, std::string fieldName, OGRFieldType fieldType);

      /** \brief Display a feature
       * @param: feature
       * @param: display geometry*/
      static void DisplayFeature(const OGRFeatureType feature, bool displayGeom = false);

      /** \brief Check if linestring is vertical or horizontal
       * @param: Linestring
       * @return true if vertical or horizontal*/
      static bool IsVerticalOrHorizontal(const OGRGeometry* geom);

      /** \brief Check if polygon is a parallelogram
       * @param :Polygone
       * @return true is parallelogram*/
      static bool IsParallelogram(const OGRPolygon* polygon);

      /** \brief Check if a geometry is valid
       * We use this method instead of geos method because a polygon can be valid even if its exterior ring is not...*/
      static bool IsValid(const OGRGeometry* geom);

      /** \brief Compute the intersection between 2 geometries and return a vector of intrsected geometries
       * @param: Geom ref
       * @param: Geom adj
       * @return :vector of geometries*/
      static std::vector<OGRGeometry*> IntersectGeoms(OGRGeometry* geomRef, OGRGeometry* geomAdj);

      /** \brief Initialize a gdal driver
       * @param: Driver name
       * @return : GDAL Driver*/
      static GDALDriver* InitializeGDALDriver(std::string driverName);

      /** \brief Clone a GDAL Dataset
       * @param: GDALDataset
       * @return :Cloned gdaldataset*/
      static GDALDataset* CloneDataset(GDALDataset* sourceDs);

      /** \brief Write OGR DS into file
       * @param: OGR Datasource
       * @param: Filename
       * @param: Layer id (-1 for all)*/
      static void WriteOGRDataSource(OGRDataSourceType* ogrDS, std::string filename, int layerId = - 1);

      /** \brief Write OGR DS into file
       * @param: OGR Datasource
       * @param: Filename
       * @param: Layer name ("" for all)*/
      static void WriteOGRDataSource(OGRDataSourceType* ogrDS, std::string filename, const std::string layerName = "");

      // \brief This method serializes a layer that can be either sent via MPI requests or written in a binary file.
      static std::vector<char> SerializeLayer(const OGRLayerType layer);

      // \brief This methods builds a layer from a bit stream. It is generic if all the needed
      // serialization methods of the specific attributes are provided by the user.
      static OGRLayerType DeSerializeLayer(const std::vector<char>& serializedLayer);


      /** \brief Convert polygon to multipolygon (used in case of self-intersecting polygon)
       * @param :Polygon
       * @return: Multipolygon (split self intersecting into multiple polygons)
       * */
      static OGRGeometry* ConvertToMultiPolygons(OGRPolygon* polygon);

      /** \brief Convert a polygon to a list of linear string (for exterior rings)
       * @param : Polygon
       * @param : Vector of linear string*/
      static std::vector<OGRLineString*> ConvertPolygonToLinearString(OGRPolygon* polygon);

      // \brief Insert a point in a linestring in order to be valid
      //@param linestring
      static void InsertPointInRightOrder(OGRLineString* linestring, OGRPoint* point);



};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaVectorOperations.txx"
#endif
