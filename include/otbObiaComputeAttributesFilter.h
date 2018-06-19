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
#ifndef otbObiaComputeAttributesFilter_h
#define otbObiaComputeAttributesFilter_h

#include <otbObiaGenericAttribute.h>
#include <otbPersistentSamplingFilterBase.h>

/**
\file otbObiaComputeAttributesFilter.h
\brief This file define the filter used to compute attributs. It extracts samples from a ROI and compute
all the attributes given
*/
namespace otb
{
namespace obia
{

/**Class specializing the merging cost function required by the generic filter*/
template <class TInputImage, class TMaskImage = otb::Image<unsigned char,2> >
class ComputeAttributesFilter : public PersistentSamplingFilterBase<TInputImage, TMaskImage>
{
	public:

		/** Standard typedefs */
	  using Self 				= ComputeAttributesFilter<TInputImage, TMaskImage>;
	  using Superclass   		= PersistentSamplingFilterBase<TInputImage, TMaskImage>;
	  using Pointer 			= itk::SmartPointer<Self>;
	  using ConstPointer	    = itk::SmartPointer<const Self>;
	  using OGRDataSourceType	= otb::ogr::DataSource;
	  using OGRDataPointer      = OGRDataSourceType::Pointer;
	  using OGRLayerType		= ogr::Layer;
	  using OGRFeatureType		= ogr::Feature;

	  /** Creation through object factory macro */
	  itkNewMacro(Self);

	  itkTypeMacro(ComputeAttributesFilter, PersistentSamplingFilterBase);


	  /** Template parameters typedefs */
	  using InputImageType 	     = TInputImage;
	  using MaskImageType 	     = TMaskImage;
	  using InputPixelType	     = typename InputImageType::InternalPixelType;
	  using RegionType 	  	     = typename InputImageType::RegionType;
	  using GenericAttributeType = GenericAttribute<InputImageType>;


	  /** Set/Get macro for the field name containing class names
	    * in the input vectors.*/
	   itkSetMacro(InputLayerName, std::string);
	   itkGetMacro(InputLayerName, std::string);

	   itkSetMacro(OutputDir, std::string);
	   itkGetMacro(OutputDir, std::string);

	   itkSetMacro(OutputFilename, std::string);
	   itkGetMacro(OutputFilename, std::string);

	   /**
	    * Reset the persistent data of the filter.
	    */
	   virtual void Reset(void);
	   /**
	    * Synthesize the persistent data of the filter.
	    */
	   virtual void Synthetize(void);

	  /**Set all computing attributs method*/
	   void SetAttributes(std::vector<GenericAttributeType*> attributes){ m_Attributes = attributes;};

	protected:
	  /** Constructor */
	  ComputeAttributesFilter();
	  /** Destructor */
	  ~ComputeAttributesFilter() ITK_OVERRIDE {};

	  /** Generate data should thread over */
	  void GenerateData(void) ITK_OVERRIDE;

	  /** Compute all attributs*/
	  void ComputeAllAttributes();

	  /**Initialize output */
	  void InitializeOutput();

	  /** Generic method called for each matching pixel position*/
	  virtual void ProcessSample(const ogr::Feature& feature,
	                             typename TInputImage::IndexType& imgIndex,
	                             typename TInputImage::PointType& imgPoint,
	                             itk::ThreadIdType& threadid) ITK_OVERRIDE;

	private:

	  //Input layer name used for computing attributs
	  std::string m_InputLayerName;

	  //Vector of GenericAttributes
	  std::vector<GenericAttributeType*> m_Attributes;

	  //Vector containing samples
	  std::map<unsigned int, std::vector<InputPixelType> > m_Samples;

	  //Input image
	  InputImageType* m_InputImage;

	  //Number of components
	  unsigned int m_NumberOfBands;

	  //Output DS
	  OGRDataSourceType::Pointer m_OutputDs;

	  //Output dir
	  std::string m_OutputDir;

	  //Output file name
	  std::string m_OutputFilename;



};

} // end of namespace obia
} // end of namespace otb
#include <otbObiaComputeAttributesFilter.txx>
#endif

