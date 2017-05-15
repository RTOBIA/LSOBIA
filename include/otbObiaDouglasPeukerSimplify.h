#ifndef otbObiaDouglasPeukerSimplify_h
#define otbObiaDouglasPeukerSimplify_h

#include "otbObiaGenericSimplifyFunc.h"

namespace otb
{
namespace obia
{

/**Class specializing the simplfy func*/
class DouglasPeukerFunc : public GenericSimplifyFunc
{
public:
    /** Some convenient alias */

    /** Standard class alias */
    using Self         = DouglasPeukerFunc;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    DouglasPeukerFunc();
    virtual ~DouglasPeukerFunc();

    /**Simplification method*/
    virtual void SimplifyLine();

    /**Set tolerance*/
    void SetTolerance(double tolerance){ m_tolerance = tolerance;};

protected:

    /**Tolerance*/
    double m_tolerance;
};

}
}
#include "otbObiaDouglasPeukerSimplify.txx"
#endif

