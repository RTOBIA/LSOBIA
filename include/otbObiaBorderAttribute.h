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
#ifndef otbObiaBorderAttribute_h
#define otbObiaBorderAttribute_h

#include "otbObiaGenericAttribute.h"

/**
\file otbObiaBorderAttribute.h
\brief This file define the Mean Attribute class used to compute attribute of all features
*/
namespace otb
{
namespace obia
{


/** \class BorderAttribute
 *	\brief Class defining the MeanAttribute. It allows to compute the mean of exterior ring with a margin.
 */
template< class TInputImage >
class BorderAttribute : public GenericAttribute<TInputImage>
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self           = BorderAttribute;
    using Pointer        = itk::SmartPointer<Self>;
    using ConstPointer   = itk::SmartPointer< const Self>;
    using PixelType	     = typename TInputImage::InternalPixelType;
    using OGRLayerType 	 = ogr::Layer;
    using OGRFeatureType = ogr::Feature;

    BorderAttribute();
    BorderAttribute(double margin);
    virtual ~BorderAttribute();

    /**Generate the required feature (for example instead of full polygon, we could only need exterior ring)*/
    virtual OGRFeatureType GenerateFeature(const OGRFeatureType inputFeature);

    /**Compute the attribut*/
    virtual double ComputeAttribut(std::vector<PixelType> samples);

    /**Set the margin*/
    void SetMargin(double margin){ m_Margin = margin;};

protected:

    /**Margin (default is 0)*/
    double m_Margin;

};

}
}
#include <otbObiaBorderAttribute.txx>
#endif
