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
#ifndef otbObiaSmallRegionsMergingFilter_h
#define otbObiaSmallRegionsMergingFilter_h
#include "itkImageRegionConstIterator.h"
#include <limits>

#include "otbObiaGraphToGraphFilter.h"

/**
\file otbObiaSmallRegionsMergingFilter.h
\brief This file define all classes required for  small regions merging (merging cost, heuristic, update attributes)
*/
namespace otb
{
namespace obia
{


/**Class specializing the merging cost function required by the generic filter*/
template< typename TCost, typename TGraph >
class SRMMergingCost
{

public:

    /** Some convenient alias */
    using ValueType        = TCost;
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = SRMMergingCost;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    //itkNewMacro(Self);
    SRMMergingCost();
    ~SRMMergingCost();

    static ValueType Max(){return std::numeric_limits<TCost>::max();};

    /** Get/Set methods */
    void SetMinimalSurface(uint32_t minimalSurface){m_MinimalSurface = minimalSurface;};

    const ValueType GetMax(){return m_MaxCost;};
    bool ComputeMergingCostsForThisNode(NodeType* curNode);
    bool ComputeMergingCostsForThisAdjNode(NodeType* curNode);
    ValueType ComputeMergingCost(NodeType* NodeIn, NodeType* NodeOut);

protected:

    /**Maximal cost of merging*/
    ValueType m_MaxCost;

    /**Minimal surface*/
    uint32_t m_MinimalSurface;

};

/**Class specializing the heuristic function required by the generic filter*/
template<typename TGraph >
class SRMHeuristic
{

public:
    /** Some convenient alias */
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = SRMHeuristic;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;


    /** Method for creation through the object factory. */
    //itkNewMacro(Self);
    SRMHeuristic();
    ~SRMHeuristic();

    NodeType* GetBestAdjacentNode(NodeType* curNode);

    /** Get/Set methods */
    void SetGraph(GraphPointerType graph){ m_Graph = graph;};
    void SetMinimalSurface(uint32_t minimalSurface){m_MinimalSurface = minimalSurface;};

protected:


    /**Compute spectral distance*/
    float ComputeSpectralDistance(NodeType* node1, NodeType* node2);

	NodeType* GetBestSmallAdjacentNode(NodeType* nodeIn);

	NodeType* GetMostSimilarNode(NodeType* node);

	bool CheckMutuality(NodeType* node_1, NodeType* node_2);

	/**Pointer to the graph*/
	GraphPointerType m_Graph;

    /** Threshold for Baatz decision*/
    float m_Threshold;

    /**Minimal surface*/
    uint32_t m_MinimalSurface;


};

template<typename TGraph >
class SRMUpdateAttribute
{

public:
    /** Some convenient alias */
    using GraphType        = TGraph;
    using GraphPointerType = typename GraphType::Pointer;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = SRMUpdateAttribute;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    //itkNewMacro(Self);
    SRMUpdateAttribute();
    ~SRMUpdateAttribute();

    /**Update attributes of nodeIn*/
    void UpdateAttributes(NodeType * nodeIn, NodeType *  nodeOut);

protected:

};



/** \class SmallRegionMergingFilter
 *    \brief Class that builds an adjacency graph where all smalls regions are merged
 *
 */
//Define the graph type : SRM (Small Region Merging)
//using GraphType = Graph< Node< SRMNodeAttribute, SRMEdgeAttribute > >;

template< typename TGraph >
class SmallRegionsMergingFilter : public GraphToGraphFilter< TGraph, TGraph >
{

public:

    /** Some convenient alias */
    using GraphType        = TGraph;
    using NodeType         = typename GraphType::NodeType;
    using EdgeType         = typename GraphType::EdgeType;

    /** Standard class alias */
    using Self         = SmallRegionsMergingFilter;
    using Superclass   = GraphToGraphFilter<GraphType, GraphType>;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SmallRegionsMergingFilter, GraphToGraphFilter);


    itkSetMacro(MinimalSurface, uint32_t);
    itkGetMacro(MinimalSurface, uint32_t);
    itkSetMacro(MergingOver, bool);
    itkGetMacro(MergingOver, bool);
    itkSetMacro(NumberOfIterations, int32_t);
    itkGetMacro(NumberOfIterations, int32_t);

protected:

    SmallRegionsMergingFilter();

    virtual ~SmallRegionsMergingFilter();

    void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

    void GenerateData();


private:

    /**Method doing one iteration of the algorithm
     * Return false if small regions has been found and merged, else true*/
    bool DoOneIteration();

    /**Method giving the most similar node*/
    NodeType * GetMostSimilarNode(NodeType * node);

    /**Method updating the merging cost of small regions*/
    void UpdateMergingCosts();

    /**Method merging nodes*/ /**TODO : Maybe move this method in a class allowing operation between nodes*/
    void Merge(NodeType* nodeIn, NodeType * nodeOut);

    /**Method updating specific attributes when merging node (like computing the new spectral mean)*/
    void UpdateSpecificAttributes(NodeType * nodeIn, NodeType * nodeOut);

    /**Method computing the spectral distance between 2 nodes*/
    float ComputeSpectralDistance(NodeType* node1, NodeType* node2);

    /**Method to check area*/
    bool HasNoSmallRegions();

    /**Method to check if graph valid*/
    bool IsGraphValid();

	/** Method to check mutuality*/
	bool CheckMutuality(NodeType* node_1, NodeType* node_2);

	/**Surface value to consider a small region*/
	uint32_t m_MinimalSurface;

    /**Boolean which flag if all small regions merged*/
    bool m_MergingOver;

    /**Max iteration*/
    int32_t m_NumberOfIterations;

};
} // end of namespace obia
} // end of namespace otb
#include "otbObiaSmallRegionsMergingFilter.txx"
#endif
