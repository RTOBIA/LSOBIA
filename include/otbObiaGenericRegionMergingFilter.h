#ifndef __otbObiaGenericRegionMergingFilter_h
#define __otbObiaGenericRegionMergingFilter_h
#include "otbObiaGraphToGraphFilter.h"

/**
\file otbObiaGenericRegionMergingFilter.h
\brief This file define the generic filter used to compute a segmentation using Baatz & Schäpe or to merge small regions
*/

namespace otb
{

namespace obia
{

/**\class GenericRegionMergingFilter otbObiaGenericRegionMergingFilter.h
 * \brief Class defining the generic filter.
 * It is templated by:
 * - An input graph type (like Baatz graph)
 * - An output graph type (generally same as input)
 * - A merging cost func, used to describe the cost between nodes in order to decide when to merge
 * - An heuristic func use to decide which node is the best adjacent node in a graph
 * - An update attribute func used to update meta data of merged node (like area, mean, etc ...)\n*/
template< typename TInputGraph,
          typename TOutputGraph,
          typename TMergingCostFunc,
          typename THeuristic,
          typename TUpdateAttributeFunc >
class GenericRegionMergingFilter : public GraphToGraphFilter<TInputGraph, TInputGraph>
{
public:

    /** Standard convenient alias */
    using Self = GenericRegionMergingFilter;
    using Superclass = GraphToGraphFilter<TInputGraph, TInputGraph>;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Some convenient alias */
    using InputGraphType  = TInputGraph;
    using OutputGraphType = TOutputGraph;
    using NodeType          = typename InputGraphType::NodeType;
    /** 
        Some convenient alias concerning the type of the function
        which computes the merging cost.
    */
    using MergingCostFunctionType = TMergingCostFunc;
    using MergingCostValueType = typename MergingCostFunctionType::ValueType;
    using HeuristicType = THeuristic;
    using UpdateAttributeFuncType = TUpdateAttributeFunc;


    /**\brief Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(GenericRegionMergingFilter, GraphToGraphFilter);

    itkSetMacro(MaxNumberOfIterations, unsigned int);
    itkGetMacro(MergingOver, bool);
    itkGetMacro(AppliedNumberOfIterations, unsigned int);

    itkGetMacro(MergingCostFunc, MergingCostFunctionType *);

    itkGetMacro(HeuristicFunc, HeuristicType *);

    itkGetMacro(UpdateAttributeFunc, UpdateAttributeFuncType *);

    /**\brief Set the merging cost func
     * \param Merging cost function*/
    void SetMergingCostFunc(MergingCostFunctionType * mergingCost){m_MergingCostFunc = mergingCost;};

    /**\brief Set the heuristic cost func
     * \param Heuristic function*/
    void SetHeuristicFunc(HeuristicType * heuristicFunc){m_HeuristicFunc = heuristicFunc;};

    /**\brief Set the update attribute func
     * \param Update attribute func*/
    void SetUpdateAttributeFunc(UpdateAttributeFuncType * updateAttribute){m_UpdateAttributeFunc = updateAttribute;};

    /**\brief Check validity of the object, it checks if the 3 required functions are set (not nullptr)*/
    virtual void CheckValidity();
protected:

    GenericRegionMergingFilter();
    virtual ~GenericRegionMergingFilter();

    void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

    void GenerateData();

    /** \briefGeneric method applying one iteration of a region merging procedure
        Can be override by inherited classes for very specific cases.
    */
    virtual bool DoOneIteration();

    /**\brief Generic method that merges two adjacent nodes.
     * \param: Node in
     * \param: Node out which will merge with node in
    */
    virtual void Merge(NodeType* nodeIn, NodeType * nodeOut);

    /**\brief Generic method computing the merging cost between pair of adjacent
        nodes in an adjacent graph.
        Can be override by inherited classes for specific cases.
    */
    virtual void ComputeMergingCosts();


private:

    /**\brief For a region merging process, there will be always a maximum number of iterations. */
    unsigned int m_MaxNumberOfIterations;

    /**\brief Real number of iterations done. */
    unsigned int m_AppliedNumberOfIterations;

    /*\brief It will be always necessary to know if there has been merges during the last iteration. */
    bool m_MergingOver;

    /** Pointers to functions that needs to be specialized */
    MergingCostFunctionType * m_MergingCostFunc;
    HeuristicType * m_HeuristicFunc;
    UpdateAttributeFuncType * m_UpdateAttributeFunc;
};


} // end of namespace obia

} // end of namespace otb

#include "otbObiaGenericRegionMergingFilter.txx"
#endif
