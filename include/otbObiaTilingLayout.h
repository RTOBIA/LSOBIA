/*
 * otbObiaTilingLayout.h
 *
 *  Created on: 21 nov. 2017
 *      Author: cresson
 */

#ifndef MODULES_REMOTE_LSOBIA_INCLUDE_OTBOBIATILINGLAYOUT_H_
#define MODULES_REMOTE_LSOBIA_INCLUDE_OTBOBIATILINGLAYOUT_H_


/**
\file otbObiaTilingLayout.h
\brief Computes the tiling layout for the first partial iterations
 */
namespace otb
{
namespace obia
{

/*
 * This function computes the maximum stability margin from the given tile size
 */
void ComputeMaximumStabilityMargin(
    /** INPUTS */
    const uint32_t width,
    const uint32_t height,

    /** OUTPUTS */
    uint32_t & niter,
    uint32_t & margin)
{
  //  itkDebugMacro(<< "Computing maximum stability margin");

  // Compute the stability margin. The naive strategy consider a margin value and a stable size equal.
  niter = 1;
  uint32_t maxMargin = std::min(width, height)/2;
  uint32_t currMargin = static_cast<uint32_t>(pow(2, niter + 1) - 2);
  margin = currMargin;

  while(currMargin < maxMargin)
    {
    margin = currMargin;
    niter++;
    currMargin = static_cast<uint32_t>(pow(2, niter + 1) - 2);
    }
  niter--;

  //  itkDebugMacro(<< "Number of iterations=" << niter << " margin=" << margin);

}


/*
 * This function computes the following parameters suited to the first partial iterations step:
 * * number of tiles in X
 * * number of tiles in Y
 * * with of regular tiles
 * * height of regular tiles
 * * the number of first partial iterations
 * * the margin for first partial iterations
 */
void ComputeTilingLayoutForFirstPartialSegmentation(
    /** INPUTS */
    const uint64_t oneNodeMemorySize,
    const uint64_t availableMemorySize,
    const uint32_t imageWidth,
    const uint32_t imageHeight,

    /** OUTPUTS */
    uint32_t & m_NbTilesX,
    uint32_t & m_NbTilesY,
    uint32_t & m_TileWidth,
    uint32_t & m_TileHeight,
    uint32_t & m_NumberOfFirstIterations,
    uint32_t & m_Margin
)
{
  // Compute the maximum number of nodes that can fit the memory
  // during first partial segmentation
  const uint64_t maximumNumberOfNodesInMemory = availableMemorySize / oneNodeMemorySize;

  // Number of nodes in the entire image
  const uint64_t nbOfNodesInImage = imageWidth * imageHeight;

  // Default layout: 1x1
  m_NbTilesX = 1;
  m_NbTilesY = 1;

  // Without margins, the number of tiles maximizing memory use
  // is equal to: nbOfNodesInImage / maximumNumberOfNodesInMemory.
  // Actually, there is tile margins. And the best scenario is to have
  // square tiles with margin = width/2, that is tiles 4x larger.
  // Hence the number of tiles maximizing memory use is 4x larger.
  uint32_t minimumNumberOfTiles = std::ceil(4.0 * ((float) nbOfNodesInImage) / ((float) maximumNumberOfNodesInMemory));

  // In the following steps, we will optimize tiling layout, starting from a number
  // of tiles equal to "minimumNumberOfTiles", up to a number of tiles equal to
  // 4 times the number of tiles (that is double rows/cols)
  uint32_t maximumNumberOfTiles = minimumNumberOfTiles * 4;

  // Search for layout which minimizes the criterion
  // The criterion is the ratio between compactness and memory usage
  // (i.e. tileWidth * tileHeight / maximumNumberOfNodesInMemory)
  float lowestCriterionValue = itk::NumericTraits<float>::max();
  for (uint32_t nbOfTiles = minimumNumberOfTiles ; nbOfTiles <= maximumNumberOfTiles ; nbOfTiles++)
    {
    // Get the multiples of k. For each one, compute the criterion of the tiling
    for (uint32_t layoutNCol = 1; layoutNCol<=nbOfTiles; layoutNCol++)
      {
#ifdef OTB_USE_MPI
      // We want number of tiles which is a multiple of the number of MPI processes
      if (nbOfTiles % layoutNCol == 0 && // Is it a multiple of the nb of Tiles and nProcs?
          nbOfTiles % otb::MPIConfig::Instance()->GetNbProcs() == 0)
#else
        if (nbOfTiles % layoutNCol == 0) // Is it a multiple of the nb of Tiles?
#endif
          {
          // Tiling layout
          uint32_t layoutNRow = nbOfTiles / layoutNCol;
          uint32_t tileWidth = imageWidth / layoutNCol;
          uint32_t tileHeight = imageHeight / layoutNRow;

          // Compute margin for regular tiles of this layout
          uint32_t maxMargin, maxIter;
          ComputeMaximumStabilityMargin(tileWidth, tileHeight, maxIter, maxMargin);
          tileWidth += 2*maxMargin;
          tileHeight += 2*maxMargin;

          // Memory use efficiency
          float percentMemory = tileWidth * tileHeight / (float) maximumNumberOfNodesInMemory; // is > 0. Could be greater than 1 in some cases!

          // Compactness
          float perimeter = tileWidth + tileHeight;
          float surface = tileWidth * tileHeight;
          float compactness = perimeter / surface * (float) std::max(tileWidth,tileHeight); // [1,+inf]

          // Update minimum criterion
          float criterion = compactness / percentMemory; // ]0, +inf]

          //          itkDebugMacro(//<< std::setprecision (2) << std::fixed
          //              << "Nb. tiles=" << nbOfTiles
          //              << " Layout: " << layoutNRow << "x" << layoutNCol
          //              << " Mem. use=" << percentMemory
          //              << " Compactness=" << compactness
          //              << " Criterion=" << criterion
          //              << " Size (no margin): " << (tileWidth-2*maxMargin)<< "x"<< (tileHeight-2*maxMargin)
          //              << " Size (with margin): " << tileWidth << "x" << tileHeight
          //              << " (margin=" << maxMargin << "/nb. iter=" << maxIter << ")" );

          if (criterion < lowestCriterionValue && percentMemory <= 1.0)
            {
            lowestCriterionValue = criterion;
            m_NbTilesX = layoutNCol;
            m_NbTilesY = layoutNRow;
            }
          }
      } // for each multiple of k
    }

  // Compute the tile size
  m_TileWidth = static_cast<uint32_t>(imageWidth/m_NbTilesX);
  m_TileHeight = static_cast<uint32_t>(imageHeight/m_NbTilesY);
  //  itkDebugMacro(<<"Selected layout: " << m_NbTilesX << "x" << m_NbTilesY
  //      << " (criterion=" << lowestCriterionValue << ")");

  // Compute the stability margin
  ComputeMaximumStabilityMargin(m_TileWidth, m_TileHeight,m_NumberOfFirstIterations, m_Margin);

  // This is the actual memory that will be used during the first partial segmentation
  uint64_t memoryUsed = availableMemorySize;
  memoryUsed *= static_cast<uint64_t>(m_TileHeight + 2*m_Margin);
  memoryUsed *= static_cast<uint64_t>(m_TileWidth + 2*m_Margin);
//  itkDebugMacro(<< "An amount of " << memoryUsed/(1024.0*1024.0) << " Mbytes of RAM will be used for regular tiles of size "
//      << (m_TileWidth + 2*m_Margin) << "x" << (m_TileHeight + 2*m_Margin) );
}

} // end namespace obia

} // end namespace otb

#endif /* MODULES_REMOTE_LSOBIA_INCLUDE_OTBOBIATILINGLAYOUT_H_ */
