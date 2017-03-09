message("OTB_MATH_PARSER MLIBS : ${OTBMathParser_LIBRARIES}")

OTB_CREATE_APPLICATION(NAME LSSegmentation
	SOURCES otbObiaLSSegmentation.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}
)
