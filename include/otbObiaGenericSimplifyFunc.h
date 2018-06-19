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
#ifndef otbObiaGenericSimplifyFunc_h
#define otbObiaGenericSimplifyFunc_h

#include <ogr_geometry.h>
//#include "itkSmartPointer.h"

namespace otb
{
namespace obia
{

/**Class specializing the simplfy func*/
class GenericSimplifyFunc
{
public:
    /** Some convenient alias */

    /** Standard class alias */
   /* using Self         = GenericSimplifyFunc;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;
*/
    /** Method for creation through the object factory. */
    GenericSimplifyFunc();
    virtual ~GenericSimplifyFunc();

    /**Simplification method*/
    virtual void SimplifyLine() = 0;

    /**Set input geometry*/
    void SetInputGeom(OGRGeometry* inputGeom){m_inputGeom = inputGeom;};

    /**Get input geometry*/
    const OGRGeometry* GetInputGeom(){return m_inputGeom;};

    /**Get outout geometry*/
    OGRGeometry* GetOutputGeom(){return m_outputGeom;};

protected:

    /**Maybe add more parameters like simplification tolerance, etc ...
     * More generaly, maybe a map containing key/value for simplification algorithm*/

    /**Input geometry*/
    OGRGeometry* m_inputGeom;

    /**Output Geometry*/
    OGRGeometry* m_outputGeom;
};

}
}

#include "otbObiaGenericSimplifyFunc.txx"
#endif

