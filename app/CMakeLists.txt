OTB_CREATE_APPLICATION(NAME LSSegmentation
	SOURCES otbObiaLSSegmentation.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} OTBCommon OTBMPIConfig
)
