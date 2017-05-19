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
      using OTBFeatureType			 = otb::ogr::Feature;
      using OTBLayerType			 = otb::ogr::Layer;

      // This method extract feature from the OGRDatasource with the field value given
      //
      // @param ogrDS: OGRDataSource
      // @param layerName : layer where the feature is
      // @param field name : name of the field to look at
      // @param field value : value of the field to look at
      // @return feature: return the feature (null if not found)
      template<typename T>
      static OTBFeatureType GetFeatureAt(OGRDataSourceType* ogrDS,
    		  	  	  	  	   	   	   	  const std::string layerName,
										  const std::string fieldName,
										  const T fieldValue);

      // This method extract all features which satisfy the field value given
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param field name: name of the field to look at
      // @param field value : value of the field to look at
      template<typename T>
      static std::vector<OTBFeatureType> GetAllFeaturesAt(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  const std::string layerName,
														  const std::string fieldName,
														  const T fieldValue);

      // This method extract all adjacent feature to the given feature
      // @param ogrDS : OGRDataSource
      // @param layerName:  layer where the features are
      // @param feature : reference feature
      // @return list of adjacent features
      static std::vector<OTBFeatureType> GetAdjacentFeatures(OGRDataSourceType* ogrDS,
	  	  	  	   	   	   	   	   	   	   	   	   	   	  	 const std::string layerName,
															 const OTBFeatureType feature);

      /**Check if a feature is valid (meaning already used for interscting for example)
       * @param: input feature
       * @return : true or false*/
      static bool isFeatureValid(const OTBFeatureType feature);

      /**Sort linestring in order to reconstruct a polygon
       * @param : Unsorted linestring
       * @return : vector of sorted linestring*/
      static std::vector<OGRLineString*> SortLinesString(std::vector<OGRLineString*> unsortedGeoms);

      // This method serializes a layer that can be either sent via MPI requests or written in a binary file.
      static std::vector<char> SerializeLayer(const OTBLayerType layer);

      // This methods builds a layer from a bit stream. It is generic if all the needed
      // serialization methods of the specific attributes are provided by the user.
      static OTBLayerType DeSerializeLayer(const std::vector<char>& serializedLayer);


};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaVectorOperations.txx"
#endif
