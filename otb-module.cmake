set(DOCUMENTATION "Large Scale OTB Application for image segmentation of very high resolution satellite scenes.")

# define the dependencies of the include module and the tests
otb_module(LSOBIA
  DEPENDS
    OTBCommon
    OTBApplicationEngine
    OTBConversion
    OTBMeanShift
    OTBMathParser
    OTBMPI
  TEST_DEPENDS
    OTBImageBase
    OTBImageIO
    OTBTestKernel
    OTBCommandLine
    OTBMathParser
  DESCRIPTION
    "Large Scale Object Based Image Analysis algorithms"
)
