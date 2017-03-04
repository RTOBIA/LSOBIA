set(DOCUMENTATION "Large Scale OTB Application for image segmentation of very high resolution satellite scenes.")

# define the dependencies of the include module and the tests
otb_module(LSOBIA
        DEPENDS
        OTBCommon
        OTBApplicationEngine
        OTBStreaming
        OTBExtendedFilename
        OTBImageIO
        OTBMeanShift
        OTBMPIConfig
        OTBConversion
        TEST_DEPENDS
        OTBTestKernel
        OTBCommandLine
        DESCRIPTION
        "${DOCUMENTATION}"
        )