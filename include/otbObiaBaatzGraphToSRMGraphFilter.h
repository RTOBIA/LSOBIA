#ifndef otbObiaBaatzGraphToSRMGraphFilter_h
#define otbObiaBaatzGraphToSRMGraphFilter_h
#include "itkImageRegionConstIterator.h"
#include <limits>

#include "otbObiaSmallRegionsMergingGraph.h"
#include "otbObiaGraphToGraphFilter.h"
#include "otbObiaImageToBaatzGraphFilter.h"

namespace otb
{
namespace obia
{
using InputGraphType  = Graph < Node < BaatzNodeAttribute, BaatzEdgeAttribute > >;
using OutputGraphType = Graph < Node < SRMNodeAttribute  , SRMEdgeAttribute   > >;

class BaatzToSRMGraphFilter : public GraphToGraphFilter<InputGraphType, OutputGraphType>
{
public:

    /** Some convenient alias */

    using InputGraphPointer    = typename InputGraphType::Pointer;
    using InputGraphConstPointer = typename InputGraphType::ConstPointer;
    using InputNodeType         = typename InputGraphType::NodeType;
    using OutputNodeType       = typename OutputGraphType::NodeType;

    /** Standard class alias */

    using Self       = BaatzToSRMGraphFilter;
    using Superclass   = GraphToGraphFilter<InputGraphType, OutputGraphType>;
    using Pointer    = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BaatzToSRMGraphFilter, GraphToGraphFilter);


protected:

    BaatzToSRMGraphFilter();
    ~BaatzToSRMGraphFilter();

    virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

    /**Convert the input graph*/
    void GenerateData();

    /**Convert a Baatz Node to SRM Node*/
    OutputNodeType* convertInputNode(InputNodeType* node);

private:

    BaatzToSRMGraphFilter(const Self &) ITK_DELETE_FUNCTION;
    void operator=(const Self &) ITK_DELETE_FUNCTION;

};
}
}
#include "otbObiaBaatzGraphToSRMGraphFilter.txx"
#endif
