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
#ifndef otbObiaSimplifyVectorFilter_h
#define otbObiaSimplifyVectorFilter_h

#include "itkProcessObject.h"
#include <otbOGRDataSourceWrapper.h>

namespace otb
{
namespace obia
{

/** \class GraphToGraphFilter
 *	\brief Base class for all process objects that requires one input graph
 	and that provides one output graph data.
 *
 */
template< typename TSimplifyFunc>
class ITK_EXPORT SimplifyVectorFilter : public itk::ProcessObject
{

public:

	/** Standard class alias */
	using Self         = SimplifyVectorFilter;
	using Superclass   = itk::ProcessObject;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Run-time type information (and related methods). */
  	itkNewMacro(Self);

  	/** Some convenient typedefs. */
  	using SimplifyFunc = TSimplifyFunc;

  	/** Return the name of the class. */
	itkTypeMacro(SimplifyVectorFilter, ProcessObject);

	  /** Definition of the input image */
	using OGRDataSourceType = otb::ogr::DataSource;
	using OGRDataSourcePointerType = OGRDataSourceType::Pointer;
	using OGRLayerType = otb::ogr::Layer;
	using OGRFeatureType = otb::ogr::Feature;
	using DataObjectPointerArraySizeType =  itk::ProcessObject::DataObjectPointerArraySizeType;

	/** Set/Get the input graph of this process object.  */
	using Superclass::SetInput;
	virtual void SetInput(const OGRDataSourceType *input);
	virtual const OGRDataSourceType * GetInput(void);

	/** Set the Field Name in which labels will be written. (default is "DN")
	* A field "FieldName" of type integer is created in the output memory layer.
	*/
	itkSetMacro(FieldName, std::string);
	/**
	* Return the Field name in which labels have been written.
	*/
	itkGetMacro(FieldName, std::string);

	/**
	* Get the layer name in which features have been written.
	*/
	itkSetMacro(LayerName, std::string);

	/**
	* Get the output \c ogr::DataSource which is a "memory" datasource.
	*/
	const OGRDataSourceType * GetOutput();

	/**
	 * Set the simplify function*/
	void SetSimplifyFunc(SimplifyFunc* simplifyFunc){m_simplifyFunc = simplifyFunc;};

	protected:

	SimplifyVectorFilter();
	~SimplifyVectorFilter() ITK_OVERRIDE;

	void GenerateInputRequestedRegion() ITK_OVERRIDE;

	/**\brief Generate Data method. Steps are following:
	 * - Initialize output dataset (setting layers definition , fields, etc ...)
	 * - Compute bounding box (to check if polygon is in the tile or not)
	 * - For each feature
	 * 		-# Get all adjacent features
	 * 		-# Compute intersection with the current feature and all adjacents and set fields to the current feature
	 * 	- Reconstruct each polygon from the edges list
	 * 	- Clean invalid polygon (for example self intersecting)
	 * 	- Remove polygons outside the tile
	 * 		*/
	void GenerateData() ITK_OVERRIDE;

	/** DataObject pointer */
	typedef itk::DataObject::Pointer DataObjectPointer;

	/**Create output*/
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

	using Superclass::MakeOutput;

	/**Initialize output DS*/
	void InitializeOutputDS();

	/**Create bounding feature to manage border*/
	OGRFeatureType* CreateBoundingBoxFeature(bool force = false);

	/**Create Tile polygon*/
	OGRPolygon* CreateTilePolygon();

	/**Intersect with bounding box and nodata layer*/
	std::vector<OGRGeometry*> IntersectWithBoundaries(OGRFeatureType feature);

	/**\brief Convert feature to edge. This method loop over all adjacents features.
	 * Then, it checks if current feature or adjacent feature is a parallelogram.
	 * If one of features is parallelogram, then we do not simplify the intersection.
	 * It also check if current feature is englobed or englobing. If it is the case, we do nothing, we just simply keep
	 * the unismplified polygon
	 * If none of the previous case, then compute intersection and simplify (if activated) the intersected line
	 * @param : Current feature
	 * @param: adjacent features
	 * @param: Layer to store englobed polygons, will be processed later*/
	void ConvertToEdges(OGRFeatureType feature, std::vector<OGRFeatureType> adjFeatures,
						OGRLayerType& insidePolygons);

	/**Create an edge
	 * This method create and add edge to map m_PolygonEdges. It calls the simplify method to simplify the intersected
	 * geometry. It is added twice: one for refFeature, and one for adjFeature (using startCoord as key)
	 * \param: Intersected geometry
	 * \param: Reference feature
	 * \param: Adjacent feature
	 * \param: Flag indicating if geometry will be simplified*/
	void AddEdge(OGRGeometry* intersectedGeometry,
				 OGRFeatureType refFeature, OGRFeatureType adjFeature,
				 bool isSimplify = true);

	/**\brief Create interior edge using englobing and englobed feature
	 * This method will allows to reconstruct interior ring when reconstructing polygon
	 * \param: Englobing feature
	 * \param: Englobed feature
	 * \param: Flag to activate simplification*/
	void CreateInteriorEdge(OGRFeatureType englobingFeature, OGRFeatureType englobedFeature, bool isSimplify = true);

	/**Compute interected features*/
	void IntersectFeatures(OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify = true);

	/**\brief Convert a geometry collection to edge. All geometries inside the geometry collection are converted
	 * to edge by using AddEdge method.
	 * \param: Intersected geometry
	 * \param: Ref Feature used for AddEdge
	 * \param: Adjacent feature used for addEdge*/
	void ConvertGeometryCollectionToEdge(OGRGeometry* intersectedLine, OGRFeatureType refFeature, OGRFeatureType adjFeature, bool isSimplify = true);

	/**\brief Reconstruct all polygons from lines geometry. A loop over all key inside the map m_PolygonEdges is done.
	 * All edges for a given key are extracted and rebuilt as a polygon by calling ReconstructPolygon*/
	void ReconstructAllPolygons();

	/** Reconstruct polygon. It first creates all linestring (made wtih 2 points), and then sort the vector containing
	 * these linestring. Once it is sorted, it creates a new geometry taking all linestring from the sorted vector as
	 * exterior rings.
	 * @param Coord of polygon
	 * @param : Vector for adjcent coords
	 * @return : created feature*/
	OGRPolygon* ReconstructPolygon(double startCoords, std::vector<double>& adjCoords);

	/** Create feature associated to polygon
	 * @param : OGRPolygon to add
	 * @param : layer to add the polygon*/
	void AddPolygon(OGRPolygon* reconstructedPolygon, double startCoords,
				    std::vector<double> adjCoords, OGRLayerType& ogrLayer);

	/**Create interior rings for polygon
	 *@param Coord of the polygon
	 *@return : vector of Linear Ring*/
	std::vector<OGRLinearRing*> CreateInteriorRings(double startCoords);

	/** Sort linestring*/
	std::vector<OGRLineString*> ConvertToGeometries(double startCoords, std::vector<double>& adjCoords);

	/**Create interior rings for fixed polygons
	 * @param: Polygon after repairing it
	 * @param: Start coords
	 * @param: New geometry with "holes" (removed englobed polygons)*/
	OGRGeometry* AddInteriorRings(OGRGeometry* fixedPolygon,double startCoords);

	/** Clean OGR Layer. This method will repear invalid polygon (like self-intersecting).*/
	void CleanLayer();

	/** Remove overlapping after fixing features.
	 *  Fixed polygons are inserted at the end of the dataset (because processed at the end), so we use
	 *  nbFixedFeatures to loop over fixed polygons in order to remove overlapping (fixed polygons tends to overlap)
	 * \param: Reconstructed layer
	 * \param: Number of fixed features*/
	void RemoveOverlappingFeatures(OGRLayerType& reconstructedLayer, unsigned int nbFixedFeatures);

	/** Remove border polygons*/
	void RemoveTileBorderPolygons(OGRLayerType& layer);

	private:

	SimplifyVectorFilter(const Self &);  //purposely not implemented

	void operator =(const Self&);      //purposely not implemented

	std::string m_FieldName;

	SimplifyFunc* m_simplifyFunc;

	std::string m_LayerName;

	/**BB Features*/
	OGRFeatureType* m_BBFeature;

	/**No data layer*/
	OGRLayerType m_NodataLayer;

	//Layer edge
	OGRLayerType m_EdgeLayer;

	//Feature def
	OGRFeatureDefn* m_EdgeFeatureDefn;

	//Current coord
	double m_CurrentCoords;

	//Edges map
	std::map<double, std::vector<OGRFeatureType>> m_PolygonEdges;

	//Interior edges
	std::map<double, std::vector<OGRGeometry*>>	m_InteriorEdges;

	//Keep track of englobing/englobed id
	std::map<double, std::vector<double>> m_EnglobedId;

	//Map to indicate if a feature is englobed or not
	std::map<double, bool> m_IsEnglobedMap;

	//BB Edges :
	//key is coordinate of the polygon
	//value is a vector containaing all intersection with boundaries
	std::map<double, std::vector<OGRGeometry*>> m_bbEdges;


};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaSimplifyVectorFilter.txx"
#endif
