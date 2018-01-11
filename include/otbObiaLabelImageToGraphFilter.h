#ifndef otbObiaLabelImageToGraphFilter_h
#define otbObiaLabelImageToGraphFilter_h
#include "otbObiaGraph.h"
#include "otbObiaImageToGraphFilter.h"
#include "itkImageRegionConstIterator.h"
#include "otbImage.h"

namespace otb
{
namespace obia
{

template< typename TLabelPixelType >
struct LabelNodeAttribute : GraphAttribute
{
    // Label of the segment
    TLabelPixelType m_Label;

    // List of internal pixels
    std::vector< CoordValueType > m_ListOfPixels;

    virtual uint64_t GetMemorySize() const
    { 
        return (sizeof(TLabelPixelType) + 
            sizeof(std::vector< CoordValueType >) + 
            CoordValueSize * m_ListOfPixels.size()); 
    }

    virtual uint64_t GetNumberOfBytesToSerialize() const 
    { 
      return stream_offset(m_Label) + stream_offset(m_ListOfPixels);
    }
    
    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const
    {
        // Serialize the label
        to_stream(stream,m_Label,position);
  
        // Serialize the number of pixels to write
	to_stream(stream,m_ListOfPixels,position);
    }

    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position)
    {
        // Deserialize the label value
        from_stream(stream,m_Label,position);

	// Deserialize the list of pixels
	from_stream(stream,m_ListOfPixels,position);
    }
};

template< typename TLabelPixelType >
class LabelImageToGraphFilter :  public ImageToGraphFilter< 
                                            otb::Image<TLabelPixelType, 2>,
                                            Graph< Node< LabelNodeAttribute<TLabelPixelType> > >
                                        >
{
public:

    /** Some convenient alias */
    using LabelPixelType               = TLabelPixelType;
    using InputImageType               = otb::Image<LabelPixelType, 2>;
    using ImageRegionConstIteratorType = itk::ImageRegionConstIterator<InputImageType>;
    using InputImagePointerType        = typename InputImageType::Pointer;
    using LabelNodeAttributeType       = LabelNodeAttribute<LabelPixelType>;
    using NodeType                     = Node< LabelNodeAttributeType >;
    using EdgeType                     = typename NodeType::EdgeType;
    using OutputGraphType              = Graph< NodeType >;

    /** Standard class alias */
    using Self         = LabelImageToGraphFilter;
    using Superclass   = ImageToGraphFilter<InputImageType, OutputGraphType>;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(LabelImageToGraphFilter, ImageToGraphFilter);

protected:

    LabelImageToGraphFilter();

    virtual ~LabelImageToGraphFilter();
    
    void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

    void GenerateData();

private:

    LabelImageToGraphFilter(const Self &) =delete;
    void operator=(const Self &) =delete;

    /** Initialize the output graph from the input label image */
    void InitOutput();

    /** Build the output graph: it is an iterative process similar to a segmentation process */
    void BuildOutput();

    /** One iteration: merge adjacent segments with the same label */
    bool DoOneIteration();

    /** Return a pointer to an adjacent node with the same label of node if it exists */
    NodeType * GetAdjacentNodeWithSameLabel(NodeType& node);

};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaLabelImageToGraphFilter.txx"
#endif
