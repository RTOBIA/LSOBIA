message("OTB_MATH_PARSER MLIBS : ${OTBMathParser_LIBRARIES}")

OTB_CREATE_APPLICATION(NAME LSSegmentation
	SOURCES otbObiaLSSegmentation.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}
)

OTB_CREATE_APPLICATION(NAME LSSmallRegionsMerging
	SOURCES otbObiaLSSmallRegionsMerging.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}
)

OTB_CREATE_APPLICATION(NAME Polygonize
	SOURCES otbObiaPolygonize.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}  
					 ${OTBGdalAdapters_LIBRARIES} ${OTBSampling_LIBRARIES}
)

OTB_CREATE_APPLICATION(NAME LSPolygonize
	SOURCES otbObiaLSPolygonize.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}  ${OTBGdalAdapters_LIBRARIES}
)

OTB_CREATE_APPLICATION(NAME ComputeAttributes
	SOURCES otbObiaComputeAttributes.cxx
	LINK_LIBRARIES ${${otb-module}_LIBRARIES} ${OTBMathParser_LIBRARIES}  ${OTBGdalAdapters_LIBRARIES}
)
