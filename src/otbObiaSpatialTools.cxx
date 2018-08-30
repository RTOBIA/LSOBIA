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
#include "otbObiaSpatialTools.h"

#include <math.h>

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

  TilingConfiguration
  SpatialTools
  ::TilePartitionningOptimizer(const int nbBands, const int width, const int height,
		  const int nbProcs, const int nbTilesPerProc, const int memPerProc, const int maxIter)
  {
	  std::cout << "Optimizing tiles size." << std::endl;
	  std::cout << "Image is " << width << "x" << height << " pixels with " << nbBands << " bands." << std::endl;
	  std::cout << "System has " << nbProcs << " proc with " << memPerProc << "Mb memory, each processing " << nbTilesPerProc << " tiles." << std::endl;

      // Compute memory used by a single pixel
	  float margin = 1.1;
	  float graphNodeSize = margin * (4 + 305 + (8 * nbBands) + 0.3 * 4 * 16 + 3* 8);
	  float imgNodeSize = (4 + 4 + 4) * nbBands;
	  float nodeSize = graphNodeSize + imgNodeSize;

	  // Image aspect ratio
	  float ratio = float(width) / float(height);
	  int nbTiles = nbProcs * nbTilesPerProc;

	  // Number of pixels per tile
	  float nbPixelsPerTile = float(width * height) / nbTiles;

	  // Estimate a target number of divisions in width from number of pixels per tile and aspect ratio
	  int targetNbTilesWidth = int(width / floor(sqrt(nbPixelsPerTile / ratio)));

      // Find the factor of nb procs closest to the target number of divisions in width
      float dist = std::abs(targetNbTilesWidth - 1);
      int bestFactor = 1;
	  for(int factor = 2 ; factor < nbTiles + 1 ; factor++)
	  {
		  if(nbTiles % factor == 0)
		  {
			  float newDist = std::abs(targetNbTilesWidth - factor);
			  if(newDist<dist)
			  {
				  bestFactor = factor;
				  dist = newDist;
			  }
		  }
	  }

      // Compute number of tiles in width and heigth
      int nbTilesWidth = bestFactor;
      int nbTilesHeight = nbTiles / nbTilesWidth;

      // Compute tile width and height
      uint32_t tileWidth = int(ceil(float(width) / nbTilesWidth));
      uint32_t tileHeight = int(ceil(float(height) / nbTilesHeight));

      uint32_t maxIterPossible = maxIter;
      // Look for optimal initial number of iteration (Baatz)
      for(int nbIter = 1 ; nbIter < maxIter ; nbIter++)
      {
		// Compute tile size with margins
		margin = std::pow(2, nbIter + 1) - 2;
		int tile_w_margin = int(floor(tileWidth + 2 * margin));
		int tile_h_margin = int(floor(tileHeight + 2 * margin));

		// Compute memory usage for this margin
		float used_memory_graph = tile_w_margin * tile_h_margin * graphNodeSize / std::pow(10, 6);
		float used_memory_img = tile_w_margin * tile_h_margin * imgNodeSize / std::pow(10, 6);
		float used_memory_margin = used_memory_img + used_memory_graph - (tileWidth * tileHeight * (graphNodeSize + imgNodeSize)) / std::pow(10, 6);

		if(used_memory_graph + used_memory_img > memPerProc)
		{
			maxIterPossible = nbIter - 1;
			break;
		}
      }

      // Return results
      TilingConfiguration result = {tileWidth, tileHeight, maxIterPossible};
      std::cout << "Proposed Partionning: " << nbTilesWidth << "x" << nbTilesHeight
    		  << "tiles of " << tileWidth << "x" << tileHeight
			  << " pixels. Maximum number of initial iteration is: " << maxIterPossible
			  << "." << std::endl;
      return result;
  }

} // end of namespace obia

} // end of namespace otb
