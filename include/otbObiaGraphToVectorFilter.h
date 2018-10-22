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
#ifndef otbObiaGraphToVectorFilter_h
#define otbObiaGraphToVectorFilter_h

#include "itkProcessObject.h"
#include <otbOGRDataSourceWrapper.h>
#include "otbImage.h"
namespace otb
{
namespace obia
{

/** \class GraphToGraphFilter
 *	\brief Base class for all process objects that requires one input graph
 	and that provides one output graph data.
 *
 */
template< typename TInputGraph>
class ITK_EXPORT GraphToVectorFilter : public itk::ProcessObject
{

public:

	/** Standard class alias */
	using Self         = GraphToVectorFilter;
	using Superclass   = itk::ProcessObject;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Run-time type information (and related methods). */
  	itkNewMacro(Self);

  	/** Some convenient typedefs. */
    using LabelPixelType		 = unsigned int;
    using LabelImageType 		 = otb::Image< LabelPixelType, 2 >;
  	using InputGraphType         = TInputGraph;
  	using InputGraphPointer      = typename InputGraphType::Pointer;
  	using InputGraphConstPointer = typename InputGraphType::ConstPointer;

  	using OTBFeatureType			 = otb::ogr::Feature;
  	using OGRLayerType			     = otb::ogr::Layer;

  	/** Return the name of the class. */
	itkTypeMacro(GraphToVectorFilter, ProcessObject);

	  /** Definition of the input image */
	using OGRDataSourceType = otb::ogr::DataSource;
	using OGRDataSourcePointerType = OGRDataSourceType::Pointer;
	using DataObjectPointerArraySizeType =  itk::ProcessObject::DataObjectPointerArraySizeType;


	/** Set/Get the input graph of this process object.  */
	using Superclass::SetInput;
	virtual void SetInput(const TInputGraph *input);
	virtual const TInputGraph * GetInput(void);


	/**Update output information*/
	virtual void UpdateOutputInformation();


	/** Set the Field Name in which labels will be written. (default is "DN")
	* A field "FieldName" of type integer is created in the output memory layer.
	*/
	itkSetMacro(FieldName, std::string);
	/**
	* Return the Field name in which labels have been written.
	*/
	itkGetMacro(FieldName, std::string);

	/**
	 *
	 */
	itkSetMacro(OutputDir, std::string);

	/**
	* Set the value of 8-connected neighborhood option used in \c GDALPolygonize
	*/
	itkSetMacro(Use8Connected, bool);
	/**
	* Get the value of 8-connected neighborhood option used in \c GDALPolygonize
	*/
	itkGetMacro(Use8Connected, bool);

	/**
	* Set the value of Xshift used for geo transform
	*/
	itkSetMacro(Xshift, unsigned int);

	/**
	* Set the value of Xshift used for geo transform
	*/
	itkSetMacro(Yshift, unsigned int);


	/**
	* Get the output \c ogr::DataSource which is a "memory" datasource.
	*/
	OGRDataSourceType * GetOutput();
	const OGRDataSourceType * GetOutput() const;

	virtual void GraftOutput(itk::DataObject *output);
	virtual void GraftOutput(const DataObjectIdentifierType & key, itk::DataObject *output);
	virtual void GraftNthOutput(unsigned int idx, itk::DataObject *output);

	protected:

	GraphToVectorFilter();
	~GraphToVectorFilter() ITK_OVERRIDE {}

	void GenerateInputRequestedRegion() ITK_OVERRIDE;

	/** Generate Data method*/
	void GenerateData() ITK_OVERRIDE;

	/** DataObject pointer */
	typedef itk::DataObject::Pointer DataObjectPointer;

	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;
	using Superclass::MakeOutput;

	/**Initialize all fields*/
	void InitializeAllFields(OGRLayerType& poLayer);


	//Clean Layer
	void CleanOGRLayer(OGRLayerType poLayer, OGRLayerType& poLayer_cleaned, OGRLayerType& polayer_nodata);

	//Check feature

	//Clean feature
	//OGRFeature* CleanFeature(OGRFeature* feature, OGRFeatureDefn* layerDefn);

//	//Self intersecting points
//	std::vector<OGRPoint> GetSelfIntersectingPoints(OGRPolygon* ogrPolygon);
//	//Clean Geometry
//	OGRPolygon* CleanSelfIntersectingPolygon(OGRPolygon* ogrPolygon);
//
//	//Create new Geometry
//	OGRPolygon* CreateSubGeomtry(OGRPointIterator* pIt, OGRPoint selfIntersectingPoint);
//
//	//Self intetrsecting
//	bool IsSelfIntersecting(OGRPoint p, std::vector<OGRPoint> listSelf);

	//Create feature for all polygons
	void CreateAllFeatures(OGRLayerType& poLayer);

	//Update the feature like adjacent polygons, etc ...
	void UpdateFeatureFields(OGRFeature* curFeature);

	//Add Feature for a given geometry
	void SetNewFields(OGRFeature* curFeature);

	private:

	GraphToVectorFilter(const Self &);  //purposely not implemented
	void operator =(const Self&);      //purposely not implemented

	//Default field name
	//Not used
	std::string m_FieldName;

	//Connexity
	bool m_Use8Connected;

	//Input graph
	InputGraphPointer m_Graph;

	//X shift
	unsigned int m_Xshift;

	//Y shift
	unsigned int m_Yshift;

	//LUT for correpsondance between label and node id
	std::vector<LabelPixelType> m_ReverseLut;

	//Output dir
	std::string m_OutputDir;
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraphToVectorFilter.txx"
#endif
