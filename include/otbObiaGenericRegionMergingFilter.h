#ifndef __otbObiaGenericRegionMergingFilter_h
#define __otbObiaGenericRegionMergingFilter_h
#include "otbObiaGraph.h"
#include "otbObiaGraphToGraphFilter.h"
#include "itkProgressReporter.h"

namespace otb
{

namespace obia
{

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
  using NodeType = typename InputGraphType::NodeType;
  using EdgeType = typename InputGraphType::EdgeType;

  /**
		Some convenient alias concerning the type of the function
		which computes the merging cost.
   */
  using MergingCostFunctionType = TMergingCostFunc;
  using MergingCostValueType = typename MergingCostFunctionType::ValueType;
  using HeuristicType = THeuristic;
  using UpdateAttributeFuncType = TUpdateAttributeFunc;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GraphToGraphFilter, GraphSource);

  itkSetMacro(MaxNumberOfIterations, unsigned int);
  itkGetMacro(MergingOver, bool);
  itkGetMacro(AppliedNumberOfIterations, unsigned int);

  itkGetMacro(MergingCostFunc, MergingCostFunctionType *);
  //itkSetMacro(MergingCostFunc, MergingCostFunctionType *);

  itkGetMacro(HeuristicFunc, HeuristicType *);
  //itkSetMacro(HeuristicFunc, HeuristicType *);

  itkGetMacro(UpdateAttributeFunc, UpdateAttributeFuncType *);
  //itkSetMacro(UpdateAttributeFunc, UpdateAttributeFuncType *);

  void SetMergingCostFunc(MergingCostFunctionType * mergingCost){m_MergingCostFunc = mergingCost;};
  void SetHeuristicFunc(HeuristicType * heuristicFunc){m_HeuristicFunc = heuristicFunc;};
  void SetUpdateAttributeFunc(UpdateAttributeFuncType * updateAttribute){m_UpdateAttributeFunc = updateAttribute;};

  /**Check validity of the object*/
  virtual void CheckValidity();
protected:

  GenericRegionMergingFilter();
  virtual ~GenericRegionMergingFilter();

  void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

  void GenerateData();

  /**
  		Generic method applying one iteration of a region merging procedure
  		Can be override by inherited classes for very specific cases.
   */
  virtual bool DoOneIteration();

  /**
  		Generic method computing the merging cost between pair of adjacent
  		nodes in an adjacent graph.
  		Can be override by inherited classes for specific cases.
   */
  virtual void ComputeMergingCosts();

private:

  /** For a region merging process, there will be always a maximum number of iterations. */
  unsigned int m_MaxNumberOfIterations;

  /** Real number of iterations done. */
  unsigned int m_AppliedNumberOfIterations;

  /** It will be always necessary to know if there has been merges during the last iteration. */
  bool m_MergingOver;

  /** Pointers to functions that needs to be specialized */
  MergingCostFunctionType * m_MergingCostFunc;
  HeuristicType * m_HeuristicFunc;
  UpdateAttributeFuncType * m_UpdateAttributeFunc;

  /** timings (benchmarks) */
  std::vector<float>         timingsValues;
  std::vector<std::string>   timingsLabels;
};


} // end of namespace obia

} // end of namespace otb

#include "otbObiaGenericRegionMergingFilter.txx"
#endif
