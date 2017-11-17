#ifndef otbObiaImageToBaatzGraphFilter_h
#define otbObiaImageToBaatzGraphFilter_h
#include "itkImageRegionConstIterator.h"
#include <limits>

#include "otbObiaGraph.h"
#include "otbObiaImageToGraphFilter.h"
#include "otbNoDataHelper.h"

/**
\file otbObiaImageToBaatzGraphFilter.h
\brief This file define classes and structure used to describe a Baatz graph
*/
namespace otb
{
namespace obia
{

struct BaatzEdgeAttribute : GraphAttribute
{
    // Merging cost
    float m_MergingCost;

    BaatzEdgeAttribute(){}

    BaatzEdgeAttribute(const BaatzEdgeAttribute& other) 
        : m_MergingCost(other.m_MergingCost)
    {}

    virtual uint64_t GetMemorySize() const;

    virtual uint64_t GetNumberOfBytesToSerialize() const;
    
    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const
    {
        (void) stream;
        (void) position;
    }
    
    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position)
    {
        (void) stream;
        (void) position;
    }
};

struct BaatzNodeAttribute : GraphAttribute
{
    // Average spectral values for each band
    std::vector< float > m_AvgSpec;

    // Standard deviation for each spectral band
    std::vector< float > m_StdSpec;

    // Area of the segment
    uint32_t m_Area;

    // Perimeter of the segment
    uint32_t m_Perimeter;

    // To indicate if the node has merged at the previous iteration
    bool m_HasPreviouslyMerged;

    BaatzNodeAttribute(){}

    BaatzNodeAttribute(const BaatzNodeAttribute& other)
    : m_AvgSpec(other.m_AvgSpec), m_StdSpec(other.m_StdSpec),
      m_Area(other.m_Area), m_Perimeter(other.m_Perimeter),
      m_HasPreviouslyMerged(other.m_HasPreviouslyMerged)
    {}

    virtual uint64_t GetMemorySize() const;

    virtual uint64_t GetNumberOfBytesToSerialize() const;
    
    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const;
    
    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position);
};

/** \class ImageToBaatzGraphFilter
 *    \brief Class that builds an adjacency graph for the Baatz & Shäpe
 *           segmentation.
 *
 */
template< typename TInputImage >
class ImageToBaatzGraphFilter : public ImageToGraphFilter<TInputImage, 
                                                          Graph< Node<
                                                                          BaatzNodeAttribute, 
                                                                          BaatzEdgeAttribute
                                                                      > 
                                                                  > 
                                                          >
{
public:

    /** Some convenient alias */
    using InputImageType               = TInputImage;
    using ImageRegionConstIteratorType = itk::ImageRegionConstIterator<InputImageType>;
    using NodeType                     = Node< BaatzNodeAttribute, BaatzEdgeAttribute >;
    using EdgeType                     = typename NodeType::EdgeType;
    using OutputGraphType              = Graph< NodeType >;
    using OutputGraphPointerType       = typename OutputGraphType::Pointer;

    /** Standard class alias */
    using Self         = ImageToBaatzGraphFilter;
    using Superclass   = ImageToGraphFilter<InputImageType, OutputGraphType>;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ImageToBaatzGraphFilter, ImageToGraphFilter);

protected:

    ImageToBaatzGraphFilter();

    virtual ~ImageToBaatzGraphFilter();

    void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

    void GenerateData();

};


} // end of namespace obia
} // end of namespace otb
#include "otbObiaImageToBaatzGraphFilter.txx"
#endif
