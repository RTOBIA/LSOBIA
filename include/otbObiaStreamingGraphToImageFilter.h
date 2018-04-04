/*
 * otbStreamingGraphToImageFilter.h
 *
 *  Created on: 6 nov. 2017
 *      Author: cresson
 */

#ifndef MODULES_REMOTE_LSGRM_INCLUDE_OTBSTREAMINGGRAPHTOIMAGEFILTER_H_
#define MODULES_REMOTE_LSGRM_INCLUDE_OTBSTREAMINGGRAPHTOIMAGEFILTER_H_

#include "itkImageSource.h"
#include "itkExceptionObject.h"
#include "itkImageRegion.h"
#include "otbObiaGraphToImageFilter.h"

// ITK FillHole filter
#include "itkGrayscaleFillholeImageFilter.h"

// Boost R-Tree
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

// Contour
#include "otbObiaContour.h"

// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace otb
{
namespace obia
{

template <typename TGraph, typename TLabelImage>
class StreamingGraphToImageFilter : public GraphToImageFilter<TGraph, TLabelImage>
{
public:
  /** Standard class typedefs. */
  typedef StreamingGraphToImageFilter                   Self;
  typedef GraphToImageFilter<TGraph, TLabelImage>       Superclass;
  typedef itk::SmartPointer<Self>                       Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StreamingGraphToImageFilter, GraphToImageFilter);

  /** Typedefs for processing */
  typedef typename TLabelImage::RegionType              RegionType;
  typedef typename TLabelImage::IndexType               IndexType;
  typedef typename TLabelImage::SizeType                SizeType;
  typedef typename TLabelImage::SpacingType             SpacingType;
  typedef typename TLabelImage::PointType               PointType;
  typedef typename TGraph::NodeType                     NodePointerType;

  typedef itk::GrayscaleFillholeImageFilter<TLabelImage,TLabelImage> FillholeFilterType;
  typedef bg::model::point<float, 2, bg::cs::cartesian> point;
  typedef bg::model::box<point> box;
  typedef std::pair<box, unsigned> value;

  /** Compute thr R-Tree */
  virtual void GenerateRTree();

  /** Prepare image allocation at the first call of the pipeline processing */
  virtual void GenerateOutputInformation(void);

  /** Does the real work. */
  virtual void GenerateData();

  /** Overloaded SetInput() method */
  virtual void SetInput(const TGraph *input);

private:

  bgi::rtree< value, bgi::quadratic<16> > rtree;

};

} // namespace obia
} // namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include <otbObiaStreamingGraphToImageFilter.txx>
#endif

#endif /* MODULES_REMOTE_LSGRM_INCLUDE_OTBSTREAMINGGRAPHTOIMAGEFILTER_H_ */

