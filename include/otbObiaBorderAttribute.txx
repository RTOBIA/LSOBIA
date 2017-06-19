#ifndef otbObiaBorderAttribute_txx
#define otbObiaBorderAttribute_txx

#include <otbObiaBorderAttribute.h>

namespace otb
{
namespace obia
{

template< class TInputImage >
BorderAttribute<TInputImage>
::BorderAttribute()
{
	this->m_Margin 		 = 0.0;
	this->m_FieldName 	 = "BorderMean";
	this->m_AttributName = "Border Mean Attribute";
	this->m_FieldType	 = OFTReal;
}

template< class TInputImage >
BorderAttribute<TInputImage>
::BorderAttribute(double margin)
{
	this->m_Margin 		 = margin;
	this->m_FieldName 	 = "BorderMean";
	this->m_AttributName = "Border Mean Attribute";
	this->m_FieldType	 = OFTReal;
}

template< class TInputImage >
BorderAttribute<TInputImage>
::~BorderAttribute()
{

}

/**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
template< class TInputImage >
typename BorderAttribute<TInputImage>::OGRFeatureType
BorderAttribute<TInputImage>
::GenerateFeature(const OGRFeatureType inputFeature)
{
	//Create a new feature
	OGRFeatureType requiredFeature(inputFeature);

	//Extract geom type
	OGRwkbGeometryType geomType = inputFeature.GetGeometry()->getGeometryType();
	OGRGeometry* outputGeom = nullptr;

	switch(geomType)
	{
		case wkbPolygon:
		{
			//Extract ext ring
			OGRPolygon* currentPolygon = (OGRPolygon*)(inputFeature.ogr().GetGeometryRef());
			outputGeom = currentPolygon->getExteriorRing()->clone();
			break;
		}

		case wkbMultiPolygon:
		{
			//Loop all polygon
			OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)(inputFeature.ogr().GetGeometryRef());

			//Declare a new geom collection
			outputGeom = new OGRGeometryCollection();

			for(unsigned int polyId = 0; polyId < multiPoly->getNumGeometries();++polyId)
			{
				//Maybe check if it is really a polygon...
				OGRPolygon* currentPolygon = (OGRPolygon*)(multiPoly->getGeometryRef(polyId));

				//Create ext ring with margins
				OGRGeometry* extRing = currentPolygon->getExteriorRing()->clone();

				//Add exterior rings
				((OGRGeometryCollection*)outputGeom)->addGeometryDirectly(extRing);
			}

			break;
		}

		case wkbGeometryCollection:
		{
			//Loop all polygon
			OGRGeometryCollection* multiGeom = (OGRGeometryCollection*)(inputFeature.ogr().GetGeometryRef());

			//Declare a new geom collection
			outputGeom = new OGRGeometryCollection();

			for(unsigned int geomId = 0; geomId < multiGeom->getNumGeometries();++geomId)
			{
				if(multiGeom->getGeometryRef(geomId)->getGeometryType() == wkbPolygon)
				{
					OGRPolygon* currentPolygon = (OGRPolygon*)(multiGeom->getGeometryRef(geomId));

					//Create ext ring with margins
					OGRGeometry* extRing = currentPolygon->getExteriorRing()->clone();

					//Add exterior rings
					((OGRGeometryCollection*)outputGeom)->addGeometryDirectly(extRing);

				}
			}
			break;
		}

		default:
		{
			//Wrong geom
			std::cerr <<"Cannont compute Border Mean attribute for "
					  << inputFeature.GetGeometry()->getGeometryName() << std::endl;
			std::cout << "inputFeature.GetGeometry() " << inputFeature.GetGeometry()->exportToGML() << std::endl;
			break;
		}
	}

	//Force to linestring, because the exploreGeometry method of PersistantSampleFilter can process only point
	//polygon, multipolygon or linestring
	requiredFeature.SetGeometry(OGRGeometryFactory::forceToLineString(outputGeom));


	return requiredFeature;
}

/**Compute the attribut*/
template< class TInputImage >
double
BorderAttribute<TInputImage>
::ComputeAttribut(std::vector<PixelType> samples)
{
	double meanValue = 0.0;
	//Loop each samples
	for(unsigned int k = 0; k < samples.size(); k++)
	{
		PixelType sample = samples[k];
		meanValue += static_cast<double>(sample);
	}

	if(samples.size() != 0)
	{
		meanValue = meanValue/samples.size();
	}

	return meanValue;
}

}//end obia
}//end otb

#endif
