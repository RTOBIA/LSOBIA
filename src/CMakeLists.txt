set(LSOBIA_SRC
		otbObiaContour.cxx
		otbObiaSpatialTools.cxx
		)
add_library(LSOBIA ${LSOBIA_SRC})
target_link_libraries(LSOBIA ${OTBCommon_LIBRARIES} OTBMPIConfig)
otb_module_target(LSOBIA)