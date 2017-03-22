#ifndef otbObiaGraph_h
#define otbObiaGraph_h
#include <algorithm>
#include <utility>
#include "otbObiaContour.h"
#include "itkProcessObject.h"

namespace otb
{
namespace obia
{

/** \struct GraphAttribute
 * \brief struct representing an abstract graph attribute
 * 
 * This is an abstract structure representing an attribute
 * of the graph (i.e, attribute of an edge or attribute of a
 * node). This structure contains C++ abstract methods that
 * need to be override by user attributes classes.
 */
struct GraphAttribute
{
    /** Returns the memory size in bytes needed to store this object in computer memory */
    virtual uint64_t GetMemorySize() const = 0;

    /** Returns the number of bytes to serialize. It may be different from the memory size
     *  since some attributes may not need to be serialized.
     */
    virtual uint64_t GetNumberOfBytesToSerialize() const = 0;

    /** Returns a stream of bits containing those attributes */
    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const = 0;

    /** Given a stream of bits, this method rebuilds the attributes */
    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position) = 0;
};

/** \struct DummyGraphAttribute
 * \brief When there is no need to specialize attributes for a node
 * or an outgoing edge, this dummy attribute is considered by default.
 */
struct DummyGraphAttribute : GraphAttribute
{
    virtual uint64_t GetMemorySize() const { return 0; }
    virtual void Serialize(std::vector<char>& stream, uint64_t& position) const 
    {
        AvoidCompilationWarings(stream, position);
    }
    virtual uint64_t GetNumberOfBytesToSerialize() const{ return 0; }
    virtual void DeSerialize(const std::vector<char>& stream, uint64_t& position)
    {
        AvoidCompilationWarings(stream, position);
    }
private:
    // Does nothing
    virtual void AvoidCompilationWarings(const std::vector<char>& stream, uint64_t& position) const
    {
        // Avoid compilation warnings...
        (void)stream;
        (void)position;
    }
};

/** \struct Edge
 * \brief struct representing an outgoing edge in OBIA
 *
 * This is the generic representation of an outgoing edge.
 * An outgoing edge contains the position of the adjacent
 * node in the graph and an object attribute that needs to
 * be specialized by the user with respect to the obia
 * application.
 */
template< typename TNode >
struct Edge
{
    // Specific attribute type of the edge is retrieved
    // from the class TNode.
    using NodeType = TNode;
    using EdgeAttributeType = typename NodeType::EdgeAttributeType;

    Edge(){}

    Edge(const Edge& other) : 
        m_TargetId(other.m_TargetId), m_Boundary(other.m_Boundary), m_Attributes(other.m_Attributes)
    {}

    uint64_t GetMemorySize() const 
    { 
        return (IdSize + UInt32Size + m_Attributes.GetMemorySize());
    }

    uint64_t GetNumberOfBytesToSerialize() const;

    // The position of the edge in the adjacency graph
    IdType m_TargetId;

    // The boundary length between both adjacent nodes
    uint32_t m_Boundary;

    // Specific attributes attached to the edge
    EdgeAttributeType m_Attributes;
};


/** \struct Node
 * \brief struct representing a node in an adjacent graph
 * in OBIA.
 *
 * This is the generic representation of an object of interest.
 * A node is identified by its position in the graph. It contains
 * the bounding box of the region it represents in the image. It
 * contains the contour of this region represented by a freeman
 * chain code. It contains a list of outgoing edges targeting
 * its spatial adjacent nodes. Finally, it contains an object
 * representing specific attributes with respect to the obia
 * process. This object needs to be specialized by the user.
 */
template< typename TNodeAttribute = DummyGraphAttribute, typename TEdgeAttribute = DummyGraphAttribute >
struct Node
{
    /** Useful alias */
    using NodeAttributeType = TNodeAttribute;
    using EdgeAttributeType = TEdgeAttribute;
    using Self              = Node< NodeAttributeType, EdgeAttributeType >;
    using EdgeType          = Edge< Self >;
    using EdgeListType      = std::vector< EdgeType >;
    using EdgeIteratorType  = typename EdgeListType::iterator;

    Node() : m_HasToBeRemoved(false), m_Valid(true) {}

    Node(const Node& other) :
         m_HasToBeRemoved(other.m_HasToBeRemoved), m_Valid(other.m_Valid),
         m_Id(other.m_Id), m_BoundingBox(other.m_BoundingBox),
         m_Contour(other.m_Contour), m_Edges(other.m_Edges),
         m_Attributes(other.m_Attributes)
    {}

    /** Returns the vectorized coordinates (x,y) of the first pixel composing this node. */
    inline uint64_t GetFirstPixelCoords() const { return m_Contour.GetStartingCoords(); }

    /** Change the vectorized coordinates (x,y) of the first pixel composing this node. */
    inline void SetFirstPixelCoords(const uint64_t newCoords){ m_Contour.SetStartingCoords(newCoords); }

    /** Given the position of the adjacent node, find the outgoing edge to it. */
    EdgeIteratorType FindEdge(const IdType targetId);

    /** Find the first edge that satisfies the lambda function predicate */
    template< typename LambdaFunctionType >
    EdgeIteratorType FindEdgeIf(LambdaFunctionType f)
    {
        return std::find_if(m_Edges.begin(), m_Edges.end(), f);
    }

    /** Add an edge and returns its address */
    inline EdgeType* AddEdge();

    /** Convenient method to apply arbitarily lambda function on each edge. */
    template< typename LambdaFunctionType >
    void ApplyForEachEdge(LambdaFunctionType f)
    {
        std::for_each(m_Edges.begin(), m_Edges.end(), f);
    }

    // Returns the memory size in bytes of the node.
    uint64_t GetMemorySize() const;

    // Return the number of bytes to serialize
    uint64_t GetNumberOfBytesToSerialize() const;

    // Flag indicating if the node has to be removed from the graph
    bool m_HasToBeRemoved:1;

    // Flag indicating if this node has to be considered for merging
    bool m_Valid:1;

    // Position of the node in the adjacency list
    IdType m_Id;

    // Bounding box (the itk::ImageRegion consumes too much memory)
    // 0: upper left x
    // 1: upper left y
    // 2: size x
    // 3: size y
    std::array<uint32_t, 4> m_BoundingBox;

    // Optimized representation of the contour
    Contour m_Contour;

    // A list of outgoing edges targeting adjacent nodes
    EdgeListType m_Edges;

    // Specific attributes related to the node
    NodeAttributeType m_Attributes;

    /* Debug */
    EdgeIteratorType Begin(){return m_Edges.begin();}
    EdgeIteratorType End(){return m_Edges.end();}
};

/** \class Graph
 * \brief class representing an adjacent graph in OBIA.
 *
 * This is the generic class for all set of object of interests
 * in the otb obia processing chain. An adjacent graph contains
 * a memory aligned list of nodes and each node containing a list of 
 * edges.
 * Each node represent an object in the image (i.e, a region during
 * a segmentation process) and contains a list of outgoing edges 
 * targeting its spatial adjacent objects.
 * This class inherits from the class itk::DataObject in ITK. However
 * Region negotiation is not applicable to the DataObject Graph. Therefore,
 * the abstract region negotiation methods are not specialized in this
 * class. 
 */

template< typename TNode >
class ITKCommon_EXPORT Graph : public itk::DataObject
{

public:

    /** Standard class typedefs */
    using Self = Graph;
    using Superclass = itk::DataObject;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self >;

    itkNewMacro(Self);
    itkTypeMacro(Graph, DataObject);

    /** Define the node and edge types */
    using NodeType         = TNode;
    using EdgeType         = typename NodeType::EdgeType;
    using NodeListType     = std::vector<NodeType>;
    using EdgeIteratorType = typename NodeType::EdgeIteratorType;

    /** Get/Set methods */
    itkSetMacro(ImageWidth, uint32_t);
    itkSetMacro(ImageHeight, uint32_t);
    itkSetMacro(NumberOfSpectralBands, uint32_t);
    itkGetConstMacro(ImageWidth, uint32_t);
    itkGetConstMacro(ImageHeight, uint32_t);
    itkGetConstMacro(NumberOfSpectralBands, uint32_t);
    itkSetMacro(ProjectionRef, std::string);
    itkGetConstMacro(ProjectionRef, std::string);

    /** To avoid multiple memory reallocation, it would be better to indicate
        as soon as possible the maximum number of nodes in the graph in order
        to reserve memory space.
    */
    void SetNumberOfNodes(const uint64_t numNodes);

    /** This method adds a node in the graph and returns its address */
    inline NodeType* AddNode();

    /** This methods initialize the node given its position in the image */
    void InitStartingNode(NodeType* node, const IdType id);

    /** This method simply returns the number of nodes in the graph */
    inline uint64_t GetNumberOfNodes() const;

    /** Given the position of the node in the graph, this methods simply
        returns the adress to this node.
     */
    inline NodeType* GetNodeAt(const IdType id);

    /** This methods removes all outgoing edges targeting this node */
    void RemoveEdgesToThisNode(NodeType& node);

    /** This method applies a lambda function on each node of the graph. */
    template< typename LambdaFunctionType >
    void ApplyForEachNode(LambdaFunctionType f);

    /** This method merges two adjacent nodes: nodeOut merges into nodeIn */
    void Merge(NodeType* nodeIn, NodeType* nodeOut);

    /** This methods removes the expired nodes from the graph, i.e
        those whose m_HasToBeRemoved is true.
    */
    void RemoveNodes();

    /** This methods returns the quantity of memory in bytes to store this graph */
    uint64_t GetMemorySize() const;

    /* Debug */
    void Print();

    using NodeIterator = typename NodeListType::iterator;
    inline NodeIterator Begin(){ return m_Nodes.begin(); }
    inline NodeIterator End(){return m_Nodes.end();}

    void InsertAtEnd(const Pointer& other)
    {
        m_Nodes.insert(m_Nodes.end(), other->Begin(), other->End());
    }

    void GraftGraphByMove(Self * other)
    {
        m_Nodes = std::move(other->m_Nodes);
        CopyAttributes(other);
    }

    void CopyGraph(ConstPointer other)
    {
        m_Nodes = other->m_Nodes;
        CopyAttributes(other);
    }

    void Reset()
    {
        m_Nodes.clear();
        m_Nodes.shrink_to_fit();
    }

protected:

    Graph(){}
    virtual ~Graph() {}
    Graph(const Self &) =delete;
    void operator=(const Self &) =delete;

    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    void CopyAttributes(ConstPointer other)
    {
        m_ImageWidth = other->m_ImageWidth;
        m_ImageHeight = other->m_ImageHeight;
        m_NumberOfSpectralBands = other->m_NumberOfSpectralBands;
        m_ProjectionRef = other->m_ProjectionRef;;
    }

    NodeListType m_Nodes;

    uint32_t m_ImageWidth;
    uint32_t m_ImageHeight;
    uint32_t m_NumberOfSpectralBands;
    std::string m_ProjectionRef;

private:

    void MergeEdge(NodeType* nodeIn, NodeType* nodeOut);
};

} // end of namespace obia
} // end of namespace otb

#include "otbObiaGraph.txx"
#endif
