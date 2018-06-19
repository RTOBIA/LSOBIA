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
#ifndef otbObiaDouglasPeukerSimplify_h
#define otbObiaDouglasPeukerSimplify_h

#include "otbObiaGenericSimplifyFunc.h"

/**
\file otbObiaDouglasPeukerSimplify.h
\brief This file define the simplify function : Douglas Peuker (implemented in GDAL)
*/
namespace otb
{
namespace obia
{

/**\class DouglasPeukerFunc otbObiaDouglasPeukerSimplify.h
 * \brief Class specializing the simplfy function*/
class DouglasPeukerFunc : public GenericSimplifyFunc
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self         = DouglasPeukerFunc;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /**\brief Constructor*/
    DouglasPeukerFunc();

    /**\brief Destructor*/
    virtual ~DouglasPeukerFunc();

    /**\brief Simplification method of input geometry
     * The method first check the type of geometry:
     * 	- If it is a multilinestring, then we force to be a linestring and GDALSimplifyWithPreserveTopology is called
     * 	- If it is already a linestring, or a point or multipoint, do nothing
     * 	- Else, the geometry given is not valid, return null geometry*/
    virtual void SimplifyLine();

    /**\brief Set the tolerance used by the GDAL simplify function*/
    void SetTolerance(double tolerance){ m_tolerance = tolerance;};

protected:

    /**\brief Tolerance used for the simplification method*/
    double m_tolerance;
};

}
}
#include "otbObiaDouglasPeukerSimplify.txx"
#endif

