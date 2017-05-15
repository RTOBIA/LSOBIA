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

