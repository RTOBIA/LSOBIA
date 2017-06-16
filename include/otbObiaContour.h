#ifndef OTB_OBIA_CONTOUR_H
#define OTB_OBIA_CONTOUR_H
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include <unordered_set>
#include "otbObiaConstExpr.h"
#include "otbObiaSpatialTools.h"


/**
\file otbObiaContour.h
\brief This file define a contour used in graph
*/
namespace otb
{

namespace obia
{

/**\class Contour otbObiaContour.h
/*	\brief Optimized representation of the contour of a region in an image using the freeman chain code.
**/
class Contour
{    
public:

    // There are only 4 possible moves along a pixel
    enum Move
    { 
        UP = 0,
        RIGHT, // 1
        DOWN, // 2
        LEFT // 3
    };

    /* A packet can contains at most 4 moves of 2 bits */
    using Packet = uint8_t;

    /* Position value type */
    using OffsetValueType = uint8_t;

    /* A contour is a list of packed moves */
    using PacketList = std::vector<Packet>;

    /* A collection of generated coordinates */
    using CoordsSet = std::unordered_set<CoordValueType>;

    /*
    *    Subclass ContourIterator
    **/
    class ContourIterator
    {

    public:

        using PacketIterator = PacketList::iterator;

        ContourIterator(PacketList *const mPtr, const OffsetValueType offset)
            : m_MovesPtr(mPtr), m_FinalOffset(offset), 
            m_CurOffset(0), m_RunningPacketIt(mPtr->begin())
        {
        }

        ~ContourIterator(){}

        inline bool IsAtEnd()
        {
            if(m_RunningPacketIt == m_MovesPtr->end())
            {
                return true;
            }
            else
            {
                auto nextIt = m_RunningPacketIt + 1;
                if(nextIt == m_MovesPtr->end() && m_CurOffset == m_FinalOffset)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }

        inline void GoToNext()
        {
            if(m_CurOffset >= 6)
            {
                ++m_RunningPacketIt;
                m_CurOffset = 0;
            }
            else
            {
                m_CurOffset += 2;
            }
        }

        inline Move GetMove()
        {
            uint8_t move = ((*m_RunningPacketIt & (3 << m_CurOffset)) >> m_CurOffset);
            return (Move)move;
        }    

    private:

        /* In order to loop over it */
        PacketList* m_MovesPtr;

        /* To determine when we reach the contour end */
        OffsetValueType m_FinalOffset;

        /* The current offset */
        OffsetValueType m_CurOffset;

        /* The running iterator over the packets */
        PacketIterator m_RunningPacketIt;
    };

    Contour() : m_Offset(8){}
    Contour(const Contour& other);
    ~Contour(){}

    inline void SetStartingCoords(const CoordValueType coords)
    {
        m_StartingCoords = coords;
    }

    inline CoordValueType GetStartingCoords() const
    {
        return m_StartingCoords;
    }

    inline void AllocateMoves(const uint64_t numMoves)
    {
        if(numMoves % 4)
        {
            m_Moves.reserve(numMoves + 1);
        }
        else
        {
            m_Moves.reserve(numMoves);
        }
    }

    /* Add one move */
    inline void AddMove(const Move move)
    {
        if(m_Offset > 6)
        {
            m_Offset = 2;
            m_Moves.push_back(move);
        }
        else
        {
            m_Moves.back() |= (move << m_Offset);
            m_Offset += 2;
        }
    }

    /* For one pixel, there are 4 moves starting from the right */
    inline void FirstInit()
    {
        AddMove(RIGHT);
        AddMove(DOWN);
        AddMove(LEFT);
        AddMove(UP);
    }

    inline ContourIterator Begin()
    {
        return ContourIterator(&m_Moves, m_Offset);
    }

    /* 
    *    Generate set of coordinates of border interior pixels
    *    along the contour.
    */
    void GenerateBorderPixels(CoordsSet& borderCells,
                            const uint32_t gridSizeX);

    /* Merge with another contour */
    void MergeWith(Contour& other, 
                   const uint32_t gridSizeX,
                   const uint32_t gridSizeY);

    /* Returns the required memory to store the contour */
    uint64_t GetMemorySize() const;

    /* Returns the number of bytes to serialize */
    uint64_t GetNumberOfBytesToSerialize() const;

    /* This methods serializes the contour into a vector of char */
    void Serialize(std::vector<char>& serializedContour, uint64_t& position) const;

    /* Given a binary vector, this methods reconstructs the contour */
    void DeSerialize(const std::vector<char>& serializedContour, uint64_t& position);

    void PrintSelf();

private:

    /* Compute eight neighborhood coordinates of the pixel located at coords */
    std::array<int64_t, 8> Get8Neighborhood(CoordValueType coords, const uint32_t gridSizeX, const uint32_t gridSizeY);

    /* 
    *    Create a new contour from the border pixels of
    *    both contours.
    */
    void CreateNewMergedContour(CoordsSet& borderCells,
                                const uint32_t gridSizeX,
                                const uint32_t gridSizeY);

    /* Starting coordinates are needed to describe a contour */
    CoordValueType m_StartingCoords;

    /* The position of the last move in the last packet */
    OffsetValueType m_Offset;

    /* A contour is a list of packed moves */
    PacketList m_Moves;

};

} // end of namespace obia

} // end of namespace otb

#endif
