/*
 * otbStreamingGraphToImageFilter.txx
 *
 *  Created on: 6 nov. 2017
 *      Author: cresson
 */

#ifndef MODULES_REMOTE_LSGRM_INCLUDE_OTBSTREAMINGGRAPHTOIMAGEFILTER_TXX_
#define MODULES_REMOTE_LSGRM_INCLUDE_OTBSTREAMINGGRAPHTOIMAGEFILTER_TXX_

#include <otbObiaStreamingGraphToImageFilter.h>

namespace otb
{
namespace obia
{

template <typename TGraph, typename TLabelImage>
void
StreamingGraphToImageFilter<TGraph, TLabelImage>
::GenerateRTree()
 {
  rtree.clear();
  unsigned int count = 0;
  auto graph = const_cast< TGraph * >( this->GetInput() );
  for(auto node = graph->Begin(); node != graph->End(); ++node)
    {
    // create a bounding box
      box b(
          point(
              node->m_BoundingBox[0],
              node->m_BoundingBox[1]),
          point(
              node->m_BoundingBox[0] + node->m_BoundingBox[2],
              node->m_BoundingBox[1] + node->m_BoundingBox[3]));

      // insert new value
      rtree.insert(std::make_pair(b, count));
      count++;
    } // next node
 }

template <typename TGraph, typename TLabelImage>
void
StreamingGraphToImageFilter<TGraph, TLabelImage>
::SetInput(const TGraph *input)
 {
  Superclass::SetInput(input);
  GenerateRTree();
 }

template <typename TGraph, typename TLabelImage>
void
StreamingGraphToImageFilter<TGraph, TLabelImage>
::GenerateOutputInformation()
 {

  // TODO: use appropriate graph methods for accessing origin, spacing, ...
  auto graph = const_cast< TGraph * >( this->GetInput() );

  // Output Largest Possible Region
  IndexType index;
  index.Fill(0);
  SizeType size;
  size[0] = graph->GetImageWidth();
  size[1] = graph->GetImageHeight();
  PointType origin;
  origin[0] = graph->GetOriginX();
  origin[1] = graph->GetOriginY();

  // TODO: add spacing in the graph
//  SpacingType spacing;
//  spacing[0] = graph->GetSpacingX();
//  spacing[1] = graph->GetSpacingY();

  RegionType outputRegion(index, size );

  // Set output informations
  TLabelImage * outputPtr = this->GetOutput();
  outputPtr->SetOrigin ( origin );
//  outputPtr->SetSpacing ( spacing );
  outputPtr->SetLargestPossibleRegion( outputRegion );
  outputPtr->SetProjectionRef(graph->GetProjectionRef() );
 }


template <typename TGraph, typename TLabelImage>
void
StreamingGraphToImageFilter<TGraph, TLabelImage>
::GenerateData()
 {
  // TODO: use R-tree only if RequestedRegion != LargestPossibleRegion (m_Graph image largest possible region)

  // Allocate the output buffer
  TLabelImage * outputPtr = this->GetOutput();
  RegionType outReqRegion = outputPtr->GetRequestedRegion();
  auto graph = const_cast< TGraph * >( this->GetInput() );
  outputPtr->SetBufferedRegion(outputPtr->GetRequestedRegion());
  outputPtr->Allocate();

  // Find nodes intersecting find the output requested region
  box query_box(
      point(
          outReqRegion.GetIndex(0),
          outReqRegion.GetIndex(1)),
      point(
          outReqRegion.GetIndex(0)+outReqRegion.GetSize(0),
          outReqRegion.GetIndex(1)+outReqRegion.GetSize(1)));
  std::vector<value> result_s;
  rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));

  // Retrieve the bounding box of the intersecting nodes (a kind of "Input requested region")
  box realBBox(query_box);
  for(auto& res : result_s)
    {
      boost::geometry::expand(realBBox, res.first);
    }
  IndexType index;
  index[0] = realBBox.min_corner().get<0>();
  index[1] = realBBox.min_corner().get<1>();
  SizeType size;
  size[0] = realBBox.max_corner().get<0>() - realBBox.min_corner().get<0>();
  size[1] = realBBox.max_corner().get<1>() - realBBox.min_corner().get<1>();
  RegionType inputRequestedRegion(index, size);

  // Generate the label image
  const typename TLabelImage::InternalPixelType noDataLabel = 0;
  typename TLabelImage::Pointer labelImage = TLabelImage::New();
  labelImage->SetRegions(inputRequestedRegion);
  labelImage->Allocate();
  labelImage->FillBuffer(noDataLabel);

  using LabelImageIterator = itk::ImageRegionIterator<TLabelImage>;
  LabelImageIterator it(labelImage, inputRequestedRegion);

  // Burn boundaries
  for(auto& res : result_s)
    {
      NodePointerType * node = graph->GetNodeAt(res.second);

      Contour::CoordsSet borderPixels;
      node->m_Contour.GenerateBorderPixels(borderPixels, graph->GetImageWidth());

      for (auto& pix: borderPixels)
        {
          index[0] = pix % graph->GetImageWidth();
          index[1] = pix / graph->GetImageWidth();
          labelImage->SetPixel(index, res.second+1); // 0 is the no-data value
        }
    }

  // Fill holes
  // WARNING: the nodes MUST BE in ascending order (e.g. Top --> Down, Left --> Right)
  // Else the holes cannot be properly filled because of local minimums.
  typename FillholeFilterType::Pointer fillFilter = FillholeFilterType::New();
  fillFilter->SetInput(labelImage);
  fillFilter->Update();

  // Keep just the stable region
  LabelImageIterator outIt(outputPtr, outReqRegion);
  LabelImageIterator inIt (fillFilter->GetOutput(), outReqRegion);
  for (inIt.GoToBegin(), outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt, ++inIt)
    outIt.Set(inIt.Get());
 }

} // namespace obia
} // namespace otb

#endif

