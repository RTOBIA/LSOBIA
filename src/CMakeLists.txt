set(LSOBIA_SRC
		otbObiaContour.cxx
		otbObiaSpatialTools.cxx
		)
add_library(LSOBIA ${LSOBIA_SRC})
target_link_libraries(LSOBIA  ${OTBMathParser_LIBRARIES} ${OTBCommon_LIBRARIES} OTBMPIConfig ${OTBGdalAdapters_LIBRARIES})
otb_module_target(LSOBIA)

