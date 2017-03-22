#include "otbObiaContour.h"

namespace otb{

    namespace obia{

        Contour::Contour(const Contour& other) :
            m_StartingCoords(other.m_StartingCoords), m_Offset(other.m_Offset), m_Moves(other.m_Moves)
        {
        }

        uint64_t Contour::GetMemorySize() const
        {
                return (CoordValueSize +
                        UInt8Size +
                        sizeof(PacketList) +
                        m_Moves.size() * UInt8Size);
        }

        uint64_t Contour::GetNumberOfBytesToSerialize() const
        {
            return (CoordValueSize +  
                    UInt64Size + // the starting coordinates and the number of packets
                    UInt8Size +
                    m_Moves.size() * UInt8Size);
        }

        void Contour::Serialize(std::vector<char>& serializedContour, uint64_t& position) const
        {
            // If there is not at least one move then, there is no contour (it is impossible)
            if(m_Moves.size() > 0)
            {
                // The serialized stream is composed of:
                // 1) the starting coords => UInt64Size bytes
                // 2) the offset value => UInt8Size bytes
                // 3) the number of packets => UInt64Size bytes
                // 4) the list of packets => m_Moves.size() * UInt8Size bytes0

                // Serialize the starting coordinates
                std::memcpy(&(serializedContour[position]), &m_StartingCoords, CoordValueSize);
                position += CoordValueSize;

                // Serialize the final offset value
                std::memcpy(&(serializedContour[position]), &m_Offset, UInt8Size);
                position += UInt8Size;

                // Serialize the number of packets
                const uint64_t packetSize = m_Moves.size();
                std::memcpy(&(serializedContour[position]), &packetSize, UInt64Size);
                position += UInt64Size;

                // Serialize the packet list
                std::memcpy(&(serializedContour[position]), &(m_Moves[0]), packetSize * UInt8Size);
                position += packetSize * UInt8Size;
            }
            else
            {
                std::cerr << "Empty contour to be serialized" << std::endl;
                exit(1);
            }
        }

        void Contour::DeSerialize(const std::vector<char>& serializedContour, uint64_t& position)
        {
            // Deserialize the starting coordinates
            std::memcpy(&m_StartingCoords, &(serializedContour[position]), CoordValueSize);
            position += CoordValueSize;

            // Deserialize the final offset value
            std::memcpy(&m_Offset, &(serializedContour[position]), UInt8Size);
            position += UInt8Size;

            // Deserialize the number of packets
            uint64_t numPackets = 0;
            std::memcpy(&numPackets, &(serializedContour[position]), UInt64Size);
            position += UInt64Size;

            // Allocate the number of packets
            m_Moves.assign(numPackets, Packet());
            std::memcpy(&(m_Moves[0]), &(serializedContour[position]), numPackets * UInt8Size);
            position += numPackets * UInt8Size;
        }

        void Contour::MergeWith(Contour& other, 
                                   const uint32_t gridSizeX,
                                   const uint32_t gridSizeY)
        {
            // Generate the border pixels of both contours
            CoordsSet borderCells;
            GenerateBorderPixels(borderCells, gridSizeX);
            other.GenerateBorderPixels(borderCells, gridSizeX);


            if(m_StartingCoords > other.m_StartingCoords)
            {
                m_StartingCoords = other.m_StartingCoords;
            }

            // Create a new merged contour
            m_Moves.clear();
            m_Offset = 8;
            CreateNewMergedContour(borderCells, gridSizeX, gridSizeY);
        }

        void Contour::CreateNewMergedContour(CoordsSet& borderCells,
                                            const uint32_t gridSizeX,
                                            const uint32_t gridSizeY)
        {
            // The first move is always to the right
            AddMove(RIGHT);

            // The prev move
            Move prev = RIGHT;

            // The current coordinates
            CoordValueType currIdx = m_StartingCoords;

            for(;;)
            {
                auto neighbors = SpatialTools::EightConnectivity(currIdx, gridSizeX, gridSizeY);

                if(prev == RIGHT)
                {
                    // Pixel at up right
                    if(neighbors[1] > -1 && 
                       neighbors[2] >-1 && 
                       borderCells.find(neighbors[1]) != borderCells.end() &&
                       borderCells.find(neighbors[2]) != borderCells.end())
                    {
                        AddMove(UP);
                        currIdx -= gridSizeX;
                                                ++currIdx;
                        prev = UP;
                    }
                    // Pixel at right
                    else if(neighbors[2] >-1 &&  borderCells.find(neighbors[2]) != borderCells.end())
                    {
                        AddMove(RIGHT);
                        ++currIdx;
                    }
                    else
                    {
                        AddMove(DOWN);
                        prev = DOWN;
                    }
                }
                else if(prev == DOWN)
                {
                    if(neighbors[4] > -1 && 
                       neighbors[3] >-1 && 
                       borderCells.find(neighbors[4]) != borderCells.end() &&
                       borderCells.find(neighbors[3]) != borderCells.end())
                    {
                        AddMove(RIGHT);
                        prev = RIGHT;
                        currIdx += (gridSizeX + 1);
                    }
                    else if(neighbors[4] > -1 &&
                            borderCells.find(neighbors[4]) != borderCells.end())
                    {
                        AddMove(DOWN);
                        currIdx += gridSizeX;
                    }
                    else
                    {
                        AddMove(LEFT);
                        prev = LEFT;
                    }
                }
                else if(prev == LEFT)
                {
                    if(neighbors[6] > -1 && 
                       neighbors[5] >-1 && 
                       borderCells.find(neighbors[6]) != borderCells.end() &&
                       borderCells.find(neighbors[5]) != borderCells.end())
                    {
                        AddMove(DOWN);
                        prev = DOWN;
                        currIdx += (gridSizeX - 1);
                    }
                    else if(neighbors[6] > -1 && 
                            borderCells.find(neighbors[6]) != borderCells.end())
                    {
                        AddMove(LEFT);
                        --currIdx;
                    }
                    // In this contexte we go up along the left edge of the pixel
                    else
                    {
                        AddMove(UP);
                        prev = UP;
                    }
                }
                else
                {
                    if(neighbors[0] > -1 && 
                       neighbors[7] >-1 && 
                       borderCells.find(neighbors[0]) != borderCells.end() &&
                       borderCells.find(neighbors[7]) != borderCells.end())
                    {
                        AddMove(LEFT);
                        currIdx -= gridSizeX;
                        --currIdx;
                        prev = LEFT;
                    }
                    else if(neighbors[0] > -1 &&
                            borderCells.find(neighbors[0]) != borderCells.end())
                    {
                        AddMove(UP);
                        currIdx -= gridSizeX;
                    }
                    else
                    {
                        if(currIdx == m_StartingCoords)
                        {
                            // we close the contour, it is done
                            break;
                        }
                        else
                        {
                            // In this contexte we go right along the top edge of the pixel
                            AddMove(RIGHT);
                            prev = RIGHT;
                        }
                    }
                }
            }
        }

        void Contour::PrintSelf()
        {
            auto cIt = Begin();
            while(!cIt.IsAtEnd())
            {
                std::cout << (int)cIt.GetMove() << "->";
                cIt.GoToNext();
            }
            std::cout << std::endl;
        }

        void Contour::GenerateBorderPixels(CoordsSet& borderCells, 
                                        const uint32_t gridSizeX)
        {    
            borderCells.insert(m_StartingCoords);

            if(m_Moves.size() > 1)
            {
                auto cIt = Begin();
                Move curr, prev = cIt.GetMove();
                cIt.GoToNext();
                CoordValueType idx = m_StartingCoords;

                while(!cIt.IsAtEnd())
                {
                    curr = cIt.GetMove();

                    if(curr == UP)
                    {
                        // Impossible case is prev = 2;
                        assert(prev != DOWN);

                        //*
                        //*
                        if(prev == UP)
                        {
                            idx -= gridSizeX;
                            borderCells.insert(idx);
                        }
            
                        //  *
                        // **
                        if(prev == RIGHT)
                        {
                            // Add the pixel at the right
                            ++idx;
                            borderCells.insert(idx);

                            // Add the pixel at the top right
                            idx = idx - gridSizeX;
                            borderCells.insert(idx);
                        }
                    }
                    else if(curr == RIGHT)
                    {
                        // Impossible case
                        assert(prev != LEFT);

                        // **
                        if(prev == RIGHT)
                        {
                            ++idx;
                            borderCells.insert(idx);
                        }

                        //*
                        //**
                        if (prev == DOWN)
                        {
                            // Add the pixel at the bottom
                            idx += gridSizeX;
                            borderCells.insert(idx);

                            // Add the pixel at the bottom right
                            ++idx;
                            borderCells.insert(idx);
                        }
            
                    }
                    else if(curr == DOWN)
                    {
                        // Impossible case
                        assert(prev != UP);

                        if(prev == DOWN)
                        {
                            idx += gridSizeX;
                            borderCells.insert(idx);
                        }
                        if(prev == LEFT)
                        {
                            // Add the pixel at the left
                            --idx;
                            borderCells.insert(idx);

                            // Add the pixel at the bottom left
                            idx += gridSizeX;
                            borderCells.insert(idx);
                        }
                    }
                    else
                    {
                        // Impossible case
                        assert(prev != RIGHT);

                        if(prev == UP)
                        {
                            // Add the pixel at top
                            idx -= gridSizeX;
                            borderCells.insert(idx);

                            // Add the pixel at top right
                            --idx;
                            borderCells.insert(idx);
                        }
                        if(prev == LEFT)
                        {
                            --idx;
                            borderCells.insert(idx);
                        }
                    }

                    prev = curr;
                    cIt.GoToNext();
                }
            }
        }

    } // end of namespace obia

} // end of namespace otb
