#ifndef otbObiaImageToGraphFilter_h
#define otbObiaImageToGraphFilter_h
#include "otbObiaGraphSource.h"

namespace otb
{
namespace obia
{

template< typename TInputImage, typename TOutputGraph >
class ImageToGraphFilter : public GraphSource< TOutputGraph >
{

public:

	/** Standard class alias */
	using Self         = ImageToGraphFilter;
	using Superclass   = GraphSource<TOutputGraph>;
	using Pointer      = itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer< const Self >;

	/** Run-time type information (and related methods). */
  itkTypeMacro(ImageToGraphFilter, GraphSource);

  	/** Some convenient typedefs. */
  	using InputImageType         = TInputImage;
  	using InputImagePointer      = typename InputImageType::Pointer;
  	using InputImageConstPointer = typename InputImageType::ConstPointer;
  	using InputImageRegionType   = typename InputImageType::RegionType;
  	using InputImagePixelType    = typename InputImageType::PixelType;

  /** ImageDimension constant */
	itkStaticConstMacro(InputImageDimension, unsigned int,
	              		  TInputImage::ImageDimension);

	virtual void SetInput(const InputImageType *input);
  virtual void SetInput(unsigned int, const TInputImage *image);

  const InputImageType * GetInput() const;
  const InputImageType * GetInput(unsigned int idx) const;

protected:

	ImageToGraphFilter();
	~ImageToGraphFilter();

	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

private:
  ImageToGraphFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};

}
} // end of namespace otb
#include "otbObiaImageToGraphFilter.txx"
#endif
