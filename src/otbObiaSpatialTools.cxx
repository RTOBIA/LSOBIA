#include "otbObiaSpatialTools.h"

namespace otb
{

namespace obia
{

  /**
   * Compute the 4 neighbors of the pixel at coords vecCoords (X) in
   * the following order:
   *  0
   * 3X1
   *  2
   */
  std::array<int64_t, 4> 
  SpatialTools
  ::FourConnectivity(const uint64_t vecCoords,
             const uint32_t imgWidth, 
             const uint32_t imgHeight)
  {
    const uint32_t x = vecCoords % imgWidth;
    const uint32_t y = vecCoords / imgWidth;
    std::array<int64_t, 4> neighbors;

    // UP
    neighbors[0] = ( y > 0 ) ? vecCoords - imgWidth : -1;

    // Right
    neighbors[1] = ( x < imgWidth - 1 ) ? vecCoords + 1 : -1;

    // Down
    neighbors[2] = ( y < imgHeight - 1 ) ? vecCoords + imgWidth : -1;

    // Left
    neighbors[3] = ( x > 0 ) ? vecCoords - 1 : -1;

    return neighbors;
  }

  /**
   * Compute the 8 neighbors of the pixel at coords vecCoords (X) in
   * the following order:
   * 701
   * 6X2
   * 543
   */
  std::array<int64_t, 8>
  SpatialTools
  ::EightConnectivity(const uint64_t vecCoords,
              const uint32_t imgWidth, 
              const uint32_t imgHeight)
  {
    const uint32_t x = vecCoords % imgWidth;
    const uint32_t y = vecCoords / imgWidth;
    std::array<int64_t, 8> neighbors;

    // UP
    neighbors[0] = ( y > 0) ? vecCoords - imgWidth : -1;

    // UP right
    neighbors[1] = ( y > 0 && x < imgWidth - 1) ? vecCoords - imgWidth + 1 : -1;

    // right
    neighbors[2] = ( x < imgWidth - 1 ) ? vecCoords + 1 : -1; 

    // bottom right
    neighbors[3] = (y < imgHeight - 1 && x < imgWidth - 1) ? vecCoords + 1 + imgWidth : -1;

    // bottom
    neighbors[4] = ( y < imgHeight - 1 ) ? vecCoords + imgWidth : -1;

    // bottom left
    neighbors[5] = ( y < imgHeight - 1 && x > 0 ) ? vecCoords + imgWidth - 1 : -1;

    // left
    neighbors[6] = ( x > 0 ) ? vecCoords - 1 : -1;

    // UP left
    neighbors[7] = ( x > 0 && y > 0 ) ? vecCoords - 1 - imgWidth : -1;

    return neighbors;
  }

  void 
  SpatialTools
  ::MergeBoundingBox(std::array<uint32_t,4>& bboxIn, 
             const std::array<uint32_t, 4>& bboxOut)
  {
    const uint32_t minUX = std::min(bboxIn[0], bboxOut[0]);
    const uint32_t minUY = std::min(bboxIn[1], bboxOut[1]);
    const uint32_t maxWidth = std::max(bboxIn[0] + bboxIn[2], bboxOut[0] + bboxOut[2]);
    const uint32_t maxHeight = std::max(bboxIn[1] + bboxIn[3], bboxOut[1] + bboxOut[3]);

    bboxIn[0] = minUX;
    bboxIn[1] = minUY;
    bboxIn[2] = maxWidth - minUX;
    bboxIn[3] = maxHeight - minUY;
  }

  std::array<uint32_t, 4> 
  SpatialTools
  ::GetMergedBoundingBox(const std::array<uint32_t,4>& bboxIn,
             const std::array<uint32_t, 4>& bboxOut)
  {
    const uint32_t minUX = std::min(bboxIn[0], bboxOut[0]);
    const uint32_t minUY = std::min(bboxIn[1], bboxOut[1]);
    const uint32_t maxWidth = std::max(bboxIn[0] + bboxIn[2], bboxOut[0] + bboxOut[2]);
    const uint32_t maxHeight = std::max(bboxIn[1] + bboxIn[3], bboxOut[1] + bboxOut[3]);

    std::array<uint32_t, 4> bboxRes;
    bboxRes[0] = minUX;
    bboxRes[1] = minUY;
    bboxRes[2] = maxWidth - minUX;
    bboxRes[3] = maxHeight - minUY;

    return bboxRes;
  }

  CoordValueType
  SpatialTools
  ::TransformPixelCoordsFromTileRefToImgRef(const CoordValueType pixCoords,
                                            const ProcessingTile& tile,
                                            const uint32_t imageWidth)
  {
    /** Compute the (x,y) coordinates of the tile */
    const uint32_t xTile = pixCoords % tile.m_Frame.GetSize(0);
    const uint32_t yTile = pixCoords / tile.m_Frame.GetSize(0);

    /** Transform the tile coordinates into the image referential */
    const uint32_t xImg = tile.m_Frame.GetIndex(0) + xTile;
    const uint32_t yImg = tile.m_Frame.GetIndex(1) + yTile;

    /** Vectorized the coordinates */
    return yImg * imageWidth + xImg;
  }

  std::array<uint32_t, 4>
  SpatialTools
  ::TransformBBoxCoordsFromTileRefToImgRef(const std::array<uint32_t, 4>& bbox,
                                           const ProcessingTile& tile)
  {
    std::array<uint32_t, 4> newBbox;
    newBbox[0] = tile.m_Frame.GetIndex(0) + bbox[0];
    newBbox[1] = tile.m_Frame.GetIndex(1) + bbox[1];
    newBbox[2] = bbox[2];
    newBbox[3] = bbox[3];

    return newBbox;
  }

  bool
  SpatialTools
  ::IsBboxInsideBoundaries(const std::array<uint32_t, 4>& bbox,
                           const uint32_t lowerCol,
                           const uint32_t lowerRow,
                           const uint32_t upperCol,
                           const uint32_t upperRow)
  {
    return (bbox[0] >= lowerCol &&
            bbox[1] >= lowerRow &&
            bbox[0] + bbox[2] - 1 <= upperCol &&
            bbox[1] + bbox[3] - 1 <= upperRow);
  }

  bool
  SpatialTools
  ::IsBboxStrictlyInsideBoundaries(const std::array<uint32_t, 4>& bbox,
                                   const uint32_t lowerCol,
                                   const uint32_t lowerRow,
                                   const uint32_t upperCol,
                                   const uint32_t upperRow)
  {
    return (bbox[0] > lowerCol &&
            bbox[1] > lowerRow &&
            bbox[0] + bbox[2] - 1 < upperCol &&
            bbox[1] + bbox[3] - 1 < upperRow);
  }

  bool
  SpatialTools
  ::IsBboxOutsideBoundaries(const std::array<uint32_t, 4>& bbox,
                            const uint32_t lowerCol,
                            const uint32_t lowerRow,
                            const uint32_t upperCol,
                            const uint32_t upperRow)
  {
    return (bbox[0] > upperCol ||
            bbox[1] > upperRow ||
            bbox[0] + bbox[2] - 1 < lowerCol ||
            bbox[1] + bbox[3] - 1 < lowerRow);
  }

  bool
  SpatialTools
  ::IsPixelAtTileBorder(const CoordValueType pixCoords,
                        const uint32_t tx,
                        const uint32_t ty,
                        const uint32_t lowerCol,
                        const uint32_t lowerRow,
                        const uint32_t upperCol,
                        const uint32_t upperRow,
                        const uint32_t nbTilesX,
                        const uint32_t nbTilesY,
                        const uint32_t imageWidth)
  {

    const uint32_t x = pixCoords % imageWidth;
    const uint32_t y = pixCoords / imageWidth;
    bool isBorder = ( (tx > 0 && x == lowerCol) ||
                    (tx < nbTilesX - 1 && x == upperCol) ||
                    (ty > 0 && y == lowerRow) ||
                    (ty < nbTilesY - 1 && y == upperRow));

    return isBorder;
  }

  bool
  SpatialTools
  ::IsPixelAtBoundaries(const CoordValueType pixCoords,
                        std::unordered_set<uint32_t>& rowBounds,
                        std::unordered_set<uint32_t>& colBounds,
                        const uint32_t inputLSImageWidth)
  {
    const uint32_t x = pixCoords % inputLSImageWidth;
    const uint32_t y = pixCoords / inputLSImageWidth;

    for(const auto& bounds : rowBounds)
    {
      
      // If the y coordinate of the pixel is equal to the boudary value.
      if(y == bounds)
      {
        return true;
      }

    } // end for(const auto& bounds : rowBounds)

    // Loop over the column borders.
    for(const auto& bounds : colBounds)
    {
      // If the x coordinate of the pixel is equal to the boudary value.
      if( x == bounds )
      {
        return true;
      }

    } // end for(const auto& bounds : colBounds)

    return false;

  }



  bool 
  SpatialTools
  ::AreFloatNumbersApproximatelyEqual(const float a, const float b)
  {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<float>::epsilon() );
  }

} // end of namespace obia

} // end of namespace otb
