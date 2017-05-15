#ifndef otbObiaSimplifyVectorFilter_txx
#define otbObiaSimplifyVectorFilter_txx

#include "otbObiaSimplifyVectorFilter.h"
#include "otbGdalDataTypeBridge.h"

#include "otbImage.h"

//gdal libraries
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "gdal_alg.h"
#include "gdal_utils.h"
#include "stdint.h" //needed for uintptr_t

namespace otb
{
namespace obia
{

template <class TSimplifyFunc>
SimplifyVectorFilter<TSimplifyFunc>
::SimplifyVectorFilter() : m_FieldName("DN")
{
   this->SetNumberOfRequiredInputs(1);
   this->SetNumberOfRequiredOutputs(1);

   GDALAllRegister();

   this->ProcessObject::SetNthOutput(0, this->MakeOutput(0) );
}


template <class TSimplifyFunc>
typename SimplifyVectorFilter<TSimplifyFunc>::DataObjectPointer
SimplifyVectorFilter<TSimplifyFunc>
::MakeOutput(DataObjectPointerArraySizeType itkNotUsed(idx))
{
  return static_cast< DataObjectPointer >(OGRDataSourceType::New().GetPointer());
}

template <class TSimplifyFunc>
const typename SimplifyVectorFilter<TSimplifyFunc>::OGRDataSourceType *
SimplifyVectorFilter<TSimplifyFunc>
::GetOutput()
{
  return static_cast< const OGRDataSourceType *>(
              this->ProcessObject::GetOutput(0));
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::SetInput(const OGRDataSourceType *input)
{
  this->Superclass::SetNthInput(0, const_cast<OGRDataSourceType *>(input));
}

template <class TSimplifyFunc>
const typename SimplifyVectorFilter<TSimplifyFunc>
::OGRDataSourceType *
SimplifyVectorFilter<TSimplifyFunc>
::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
    {
    return ITK_NULLPTR;
    }

  return static_cast<const OGRDataSourceType *>(this->Superclass::GetInput(0));
}

template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::GenerateInputRequestedRegion(void)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename OGRDataSourceType::Pointer input  =
    const_cast<OGRDataSourceType *> (this->GetInput());

  if ( !input )
    {
    return;
    }
  // The input is necessarily the largest possible region.
  input->SetRequestedRegionToLargestPossibleRegion();

}


template <class TSimplifyFunc>
void
SimplifyVectorFilter<TSimplifyFunc>
::GenerateData(void)
{

	//Loop accross all geometry
	std::cout << "Generate Data for SimplifyVectorFilter" << std::endl;

}

template <class TSimplifyFunc>
GDALDriver*
SimplifyVectorFilter<TSimplifyFunc>
::initializeGDALDriver(std::string driverName)
{
	GDALDriver *poDriver;
	GDALAllRegister();
	poDriver = (GDALDriver*) GDALGetDriverByName(driverName.c_str());
	if( poDriver == NULL )
	{
		printf( "%s driver not available\n", driverName);
		exit( 1 );
	}

	return poDriver;

}


template <class TSimplifyFunc>
void SimplifyVectorFilter<TSimplifyFunc>
::WriteFile(std::string filepath)
{
	//Get output
	auto ogrDS = this->GetOutput();
	otb::ogr::Layer layer = ogrDS->GetLayer(0);
	//TODO : Temporaire
	GDALDriver*  poDriverGml = initializeGDALDriver("GML");
	std::string  output_gml = "/space/USERS/isnard/tmp/mypolygons.gml";
	GDALDataset* polyGDs = poDriverGml->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, NULL );
	polyGDs->CopyLayer(&(layer.ogr()), "test_2", nullptr);

	//Clean memory
	//Maybe see if we need to delete driver...
	GDALClose(polyGDs);
}
} //End namespace obia
} // end namespace otb

#endif
