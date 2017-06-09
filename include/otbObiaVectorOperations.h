#ifndef otbObiaVectorOperations_h
#define otbObiaVectorOperations_h

#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <unordered_set>

#include <otbOGRDataSourceWrapper.h>


namespace otb
{
namespace obia
{

    class VectorOperations
    {

      public:

      /* Some convenient alias */
      using OGRDataSourceType		 = otb::ogr::DataSource;
      using OGRFeatureType			 = otb::ogr::Feature;
      using OGRLayerType			 = otb::ogr::Layer;

      // This method extract feature from the OGRDatasource with the field value given
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

      // This method extract all features which satisfy the field value given
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param field name: name of the field to look at
      // @param field value : value of the field to look at
      template<typename T>
      static std::vector<OGRFeatureType> GetAllFeaturesAt(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  const std::string layerName,
														  const std::string fieldName,
														  const T fieldValue);

      // This method extract all adjacent feature to the given feature
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param feature : reference feature
      // @return list of adjacent features
      static std::vector<OGRFeatureType> GetAdjacentFeatures(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  	 const std::string layerName,
															 const OGRFeatureType feature);

      /**Check if a feature is valid (meaning already used for interscting for example)
       * @param: input feature
       * @return : true or false*/
      static bool isFeatureValid(const OGRFeatureType feature);

      /**Sort linestring in order to reconstruct a polygon
       * @param : Unsorted linestring
       * @return : vector of sorted linestring*/
      static std::vector<OGRLineString*> SortLinesString(std::vector<OGRLineString*> unsortedGeoms);

      /**Cast geometry to polygon geometry
       * @param: Geometry
       * @return: casted polygon*/
      static OGRPolygon* CastToPolygon(OGRGeometry* geomRef);

      /**Create a new field to the input layer*/
      static void CreateNewField(OGRLayerType& poLayer, std::string fieldName, OGRFieldType fieldType);

      /**Display a feature
       * @param: feature
       * @param: display geometry*/
      static void DisplayFeature(const OGRFeatureType feature, bool displayGeom = false);

      /**Check if linestring is vertical or horizontal
       * @param: Linestring
       * @return true if vertical or horizontal*/
      static bool IsVerticalOrHorizontal(const OGRGeometry* geom);

      /**Check if polygon is a parallelogram
       * @param :Polygone
       * @return true is parallelogram*/
      static bool IsParallelogram(const OGRPolygon* polygon);

      /**Check if a geometry is valid
       * We use this method instead of geos method because a polygon can be valid even if its exterior ring is not...*/
      static bool IsValid(const OGRGeometry* geom);

      /**Compute the intersection between 2 geometries and return a vector of intrsected geometries
       * @param: Geom ref
       * @param: Geom adj
       * @return :vector of geometries*/
      static std::vector<OGRGeometry*> IntersectGeoms(OGRGeometry* geomRef, OGRGeometry* geomAdj);

      /**Initialize a gdal driver
       * @param: Driver name
       * @return : GDAL Driver*/
      static GDALDriver* InitializeGDALDriver(std::string driverName);

      /**Clone a GDAL Dataset
       * @param: GDALDataset
       * @return :Cloned gdaldataset*/
      static GDALDataset* CloneDataset(GDALDataset* sourceDs);

      /**Write OGR DS into file
       * @param: OGR Datasource
       * @param: Filename
       * @param: Layer id (-1 for all)*/
      static void WriteOGRDataSource(OGRDataSourceType* ogrDS, std::string filename, int layerId = - 1);

      // This method serializes a layer that can be either sent via MPI requests or written in a binary file.
      static std::vector<char> SerializeLayer(const OGRLayerType layer);

      // This methods builds a layer from a bit stream. It is generic if all the needed
      // serialization methods of the specific attributes are provided by the user.
      static OGRLayerType DeSerializeLayer(const std::vector<char>& serializedLayer);


      /**Convert polygon to multipolygon (used in case of self-intersecting polygon)
       * @param :Polygon
       * @return: Multipolygon (split self intersecting into multiple polygons)
       * */
      static OGRGeometry* ConvertToMultiPolygons(OGRPolygon* polygon);

      /**Add self intersecting point into a linear ring used to recreate a polygon with more points
       * @param: Input polygon
       * @return : New polygon with ring */
      static OGRPolygon* CreatePolygonWithSelfIntersectingPoints(OGRPolygon* polygon);
//
//      /**Clean feature by removing self intersecting point
//       * @param: Feature to clean
//       * @param: Layer definition
//       * @param: Field name (used to identify cleaned feeature)
//       * @return : a cleaned feature*/
//      static OGRFeature* CleanFeature(OGRFeature* feature, OGRFeatureDefn* layerDefn, std::string fieldName = otb::obia::labelFieldName);
//
//      /**Remove duplicated edges from a polygon.
//       * @param: Polygon to clean
//       * @return: Cleaned polygon*/
//      static OGRPolygon* RemoveDuplicatedEdges(OGRPolygon* polygon);
//
//
//      protected:
//
//      //Clean Geometry
//      static OGRPolygon* CleanSelfIntersectingPolygon(OGRPolygon* ogrPolygon);
//
//      //Create new Geometry
//      static OGRPolygon* CreateSubGeomtry(OGRPointIterator* pIt, OGRPoint selfIntersectingPoint);
//
//      //Self intersecting points
//      static std::vector<OGRPoint> GetSelfIntersectingPoints(OGRPolygon* ogrPolygon);
//
//      //Self intetrsecting
//      static bool IsSelfIntersecting(OGRPoint p, std::vector<OGRPoint> listSelf);
//
      /**Convert a polygon to a list of linear string (for exterior rings)
       * @param : Polygon
       * @param : Vector of linear string*/
      static std::vector<OGRLineString*> ConvertPolygonToLinearString(OGRPolygon* polygon);

      //Insert a point in a linestring in order to be valid
      //@param linestring
      static void InsertPointInRightOrder(OGRLineString* linestring, OGRPoint* point);



};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaVectorOperations.txx"
#endif
