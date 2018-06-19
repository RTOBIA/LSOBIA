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
#ifndef otbObiaGraphToImageFilter_h
#define otbObiaGraphToImageFilter_h
#include "itkImageSource.h"
#include "itkImageRegionIterator.h"

namespace otb
{
namespace obia
{

template< typename TInputGraph, typename TOutputImage >
class GraphToImageFilter : public itk::ImageSource<TOutputImage>
{

public:

    /** Standard class alias */
    using Self         = GraphToImageFilter;
    using Superclass   = itk::ImageSource<TOutputImage>;
    using Pointer      = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer< const Self >;

    /** Run-time type information (and related methods). */
    itkTypeMacro(GraphToImageFilter, ImageSource);

    /** Some convenient typedefs. */
    using InputGraphType         = TInputGraph;
    using InputGraphPointer      = typename InputGraphType::Pointer;
    using InputGraphConstPointer = typename InputGraphType::ConstPointer;
    using OutputImageType           = TOutputImage;
    using OutputImagePointerType  = typename TOutputImage::Pointer;
    using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;

    virtual void SetInput(const InputGraphType *input);
    virtual void SetInput(unsigned int, const InputGraphType *image);

    const InputGraphType * GetInput() const;
    const InputGraphType * GetInput(unsigned int idx) const;

protected:

    GraphToImageFilter();
    ~GraphToImageFilter();
    virtual void UpdateOutputInformation();
    virtual void GenerateOutputInformation();
    virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

private:
    GraphToImageFilter(const Self &) =delete;
    void operator=(const Self &) =delete;

};

} // end of namespace obia
} // end of namespace otb
#include "otbObiaGraphToImageFilter.txx"
#endif
