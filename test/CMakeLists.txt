if (BUILD_TESTING)
	otb_module_test()
	# Handle multiple build environments
	if(NOT OTB_SOURCE_DIR)
	    set(TEMP_DIR ${LSOBIA_BINARY_DIR}/Testing/Temporary/)
	else()
	    set(TEMP_DIR ${TEMP})
	endif()
endif()

# ----------------------------------------------------------------------------------------------
# TODO : ADD: memory check tests (valgrind)?

# ----------------------------------------------------------------------------------------------
# Segmentation Application Tests
# - Baatz Tests
otb_test_application(NAME apLSOBIA_BaatzOneProc
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "BaatzOneProcLabel"
	-io.temp ${TEMP_DIR}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 5
	-algorithm.baatz.numitpartial 5
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-algorithm.baatz.aggregategraphs "on"
	-processing.writeimages "on"
	-processing.writegraphs "on"
	-processing.memory 2000
	-processing.maxtilesizex 1000
	-processing.maxtilesizey 1000
	)

otb_test_application(NAME apLSOBIA_BaatzOneProcTiled
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "BaatzOneProcTiled"
	-io.temp ${TEMP_DIR}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 5
	-algorithm.baatz.numitpartial 5
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-algorithm.baatz.aggregategraphs "on"
	-processing.writeimages "on"
	-processing.writegraphs "on"
	-processing.memory 2000
	-processing.maxtilesizex 250
	-processing.maxtilesizey 250
	)

otb_test_application(NAME apLSOBIA_BaatzOneProcLittleTiled
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_100_100.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "BaatzOneProcLabel"
	-io.temp ${TEMP_DIR}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 3
	-algorithm.baatz.numitpartial 5
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-algorithm.baatz.aggregategraphs "on"
	-processing.writeimages "on"
	-processing.writegraphs "on"
	-processing.memory 200
	-processing.maxtilesizex 25
	-processing.maxtilesizey 25
	)
	
otb_test_application(NAME apLSOBIA_BaatzOneProcLittle
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_100_100.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "BaatzOneProcLabel"
	-io.temp ${TEMP_DIR}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 1
	-algorithm.baatz.numitpartial 1
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-algorithm.baatz.aggregategraphs "on"
	-processing.writeimages "on"
	-processing.writegraphs "on"
	-processing.memory 2000
	-processing.maxtilesizex 100
	-processing.maxtilesizey 100
	)
	
otb_test_application(NAME apLSOBIA_GraphPolygonizeOneProc
	APP GraphPolygonize
	OPTIONS -io.gr ${LSOBIA_SOURCE_DIR}/data/Graph_little.dat
	-io.out.dir ${TEMP_DIR}
	-io.out.gmlfile "PolygonsFile.gml"
	-io.temp ${TEMP_DIR}
	)	
	
otb_test_application(NAME apLSOBIA_LSPolygonizeOneProc
	APP LSPolygonize
	OPTIONS -io.gr ${LSOBIA_SOURCE_DIR}/data/Graph_0_0.dat
	-io.out.dir ${TEMP_DIR}
	-io.out.gmlfile "PolygonsFile.gml"
	-io.temp ${TEMP_DIR}
)	

# ----------------------------------------------------------------------------------------------
# Segmentation Application Tests
# - Mean Shift Tests
otb_test_application(NAME apLSOBIA_MeanShiftOneProc
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "MeanShiftOneProcLabel"
	-io.temp ${TEMP_DIR}
	-algorithm meanshift
	-algorithm.meanshift.maxiter 100
	-algorithm.meanshift.spatialr 15
	-algorithm.meanshift.spectralr 16
	-algorithm.meanshift.threshold 0.01
	-algorithm.meanshift.ranger 0
	-algorithm.meanshift.modesearch "off"
	-processing.memory 2000
	-processing.maxtilesizex 1000
	-processing.maxtilesizey 1000)

# ----------------------------------------------------------------------------------------------
set(${otb-module}TestDriver
	otbLSOBIATestDriver.cxx
	otbObiaStreamUtils.cxx
	otbObiaConvertGraphToImage.cxx
)	

add_executable(otb${otb-module}TestDriver ${${otb-module}TestDriver})
target_link_libraries(otb${otb-module}TestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otb${otb-module}TestDriver)

otb_add_test(NAME otbObiaStreamUtilsTest COMMAND otbLSOBIATestDriver otbObiaStreamUtils)


# ----------------------------------------------------------------------------------------------
# Graph to image conversion test driver
if(OTB_USE_GDAL_20)
	set(${otb-module}ConvertGraphToImageTests
	    otbLSOBIAConvertGraphToImageTestDriver.cxx
	    otbObiaConvertGraphToImage.cxx
	    )

	add_executable(otbLSOBIAConvertGraphToImageTestDriver ${${otb-module}ConvertGraphToImageTests})
	target_link_libraries(otbLSOBIAConvertGraphToImageTestDriver ${${otb-module}-Test_LIBRARIES})
	otb_module_target_label(otbLSOBIAConvertGraphToImageTestDriver)

	otb_add_test(NAME otbLSOBIAConvert COMMAND otbLSOBIAConvertGraphToImageTestDriver
		     otbObiaConvertGraphToImage
		     ${LSOBIA_SOURCE_DIR}/data/Graph_0_1.dat
		     ${TEMP_DIR}/image_from_graph.tif)
else(OTB_USE_GDAL_20)
  MESSAGE("OTB does not use a GDAL version > 2.0, not building GDAL tests")
endif(OTB_USE_GDAL_20)


# ----------------------------------------------------------------------------------------------
# Large Scale Baatz Segmentation test driver
set(${otb-module}LSBaatzSegmentationTests
    otbLSOBIALSBaatzSegmentationTestDriver.cxx
    otbObiaLSBaatzSegmentation.cxx
    )

add_executable(otbLSOBIALSBaatzSegmentationTestDriver ${${otb-module}LSBaatzSegmentationTests})
target_link_libraries(otbLSOBIALSBaatzSegmentationTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIALSBaatzSegmentationTestDriver)

otb_add_test(NAME otbLSOBIALSBaatzSegmentation COMMAND otbLSOBIALSBaatzSegmentationTestDriver
             otbObiaLSBaatzSegmentation
             ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
             1000
             1000
             5
             5
             ${TEMP_DIR}
             2000
             40
             0.05
             0.95
             ${TEMP_DIR}
	     	 "on"
             "on"
             "on"
             "LSBaatzSegmentationTestLabelImage"
)


otb_add_test_mpi(NAME otbLSOBIALSBaatzSegmentationMPI COMMAND otbLSOBIALSBaatzSegmentationTestDriver
             NBPROCS 4
             otbObiaLSBaatzSegmentation
             ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
             250
             250
             5
             5
             ${TEMP_DIR}
             2000
             40
             0.05
             0.95
             ${TEMP_DIR}
	         "on"
             "on"
             "on"
             "LSBaatzSegmentationMPITestLabelImage"
	     )


# ----------------------------------------------------------------------------------------------
# Large Scale Mean Shift test driver
set(${otb-module}LSMeanShiftTests
    otbLSOBIALSMeanShiftTestDriver.cxx
    otbObiaLSMeanShift.cxx
    )

add_executable(otbLSOBIALSMeanShiftTestDriver ${${otb-module}LSMeanShiftTests})
target_link_libraries(otbLSOBIALSMeanShiftTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIALSMeanShiftTestDriver)

otb_add_test(NAME otbLSOBIALSMeanShift COMMAND otbLSOBIALSMeanShiftTestDriver
             otbObiaLSMeanShift
             ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
             1000
             1000
             2000
             ${TEMP_DIR}
             100
             15
             16
             0.1
             0
             0
             ${TEMP_DIR}
             "LSMeanShiftLabelImage"
	     )

otb_add_test_mpi(NAME otbLSOBIALSMeanShiftMPI COMMAND otbLSOBIALSMeanShiftTestDriver
             NBPROCS 4
             otbObiaLSMeanShift
             ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
             1000
             1000
             2000
             ${TEMP_DIR}
             100
             15
             16
             0.1
             0
             0
             ${TEMP_DIR}
             "LSMeanShiftMPILabelImage"
	     )


# ----------------------------------------------------------------------------------------------
# Large Scale Small Regions Merging test driver
# TODO : add mpi test
set(${otb-module}LSSmallRegionsMergingTests
    otbLSOBIALSSmallRegionsMergingTestDriver.cxx
    otbObiaLSSmallRegionsMerging.cxx
    )

add_executable(otbLSOBIALSSmallRegionsMergingTestDriver ${${otb-module}LSSmallRegionsMergingTests})
target_link_libraries(otbLSOBIALSSmallRegionsMergingTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIALSSmallRegionsMergingTestDriver)

# ----------------------------------------------------------------------------------------------
# Graph reader test driver (same functions as Graph to image conversion test driver)
set(${otb-module}ReadGraphTests
    otbLSOBIAReadGraphTestDriver.cxx
    otbObiaReadGraph.cxx
    )

add_executable(otbLSOBIAReadGraphTestDriver ${${otb-module}ReadGraphTests})
target_link_libraries(otbLSOBIAReadGraphTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIAReadGraphTestDriver)

otb_add_test(NAME otbLSOBIAReadGraph COMMAND otbLSOBIAReadGraphTestDriver
             otbObiaReadGraph
             ${LSOBIA_SOURCE_DIR}/data/Graph_0_1.dat
	     ${TEMP_DIR}/image_from_graph_without_gdal.tif)


# ----------------------------------------------------------------------------------------------
# Small Region Merging on a single file test driver
set(${otb-module}SmallRegionMergingOnFileTests
    otbLSOBIASmallRegionMergingOnFileTestDriver.cxx
    otbObiaSmallRegionMergingOnFile.cxx
    )

add_executable(otbLSOBIASmallRegionMergingOnFileTestDriver ${${otb-module}SmallRegionMergingOnFileTests})
target_link_libraries(otbLSOBIASmallRegionMergingOnFileTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIASmallRegionMergingOnFileTestDriver)

# TODO : add test when SmallRegionMerging is operational
#otb_add_test(NAME otbLSOBIASmallRegionMergingOnFile COMMAND otbLSOBIASmallRegionMergingOnFileTestDriver
#             otbObiaSmallRegionMergingOnFile
#             ${LSOBIA_SOURCE_DIR}/data/Graph_0_0.dat
#	     500
#             ${TEMP_DIR}
#             )


# ----------------------------------------------------------------------------------------------
# Graph Splitting test driver
set(${otb-module}SplitGraphTests
    otbLSOBIASplitGraphTestDriver.cxx
    otbObiaSplitGraph.cxx
    )

add_executable(otbLSOBIASplitGraphTestDriver ${${otb-module}SplitGraphTests})
target_link_libraries(otbLSOBIASplitGraphTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIASplitGraphTestDriver)

otb_add_test(NAME otbLSOBIASplitGraph COMMAND otbLSOBIASplitGraphTestDriver
             otbObiaSplitGraph
             ${LSOBIA_SOURCE_DIR}/data/Graph_0_0.dat
             ${TEMP_DIR}
             250
             250
             )

# ----------------------------------------------------------------------------------------------
# Small Regions Merging Test

otb_test_application(NAME apLSOBIA_SRMOneProc
	APP LSSmallRegionsMerging
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out.dir ${TEMP_DIR}
	-io.out.labelimage "SRMOneProcLabelImage"
	-io.temp ${TEMP_DIR}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 5
	-algorithm.baatz.numitpartial 5
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-processing.memory 2000
	-processing.maxtilesizex 500
	-processing.maxtilesizey 500
	-processing.minimalsurface 200
	)

#---------------------------------------------------------------------------------------------------
# Polygonize Test

# Large Scale Polygonize test driver
set(${otb-module}LSPolygonizeTests
    otbLSOBIALSPolygonizeTestDriver.cxx
    otbObiaLSPolygonizeTest.cxx
    )

add_executable(otbLSOBIALSPolygonizeTestDriver ${${otb-module}LSPolygonizeTests})
target_link_libraries(otbLSOBIALSPolygonizeTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIALSPolygonizeTestDriver)

otb_add_test_mpi(NAME otbObiaLSPolygonizeMPI COMMAND otbLSOBIALSPolygonizeTestDriver
             NBPROCS 4
             otbObiaLSPolygonize
             ${LSOBIA_SOURCE_DIR}/data/
             ${TEMP_DIR}
             ${TEMP_DIR}
             1000
             1000
             2
             2
             500
             500
	     )
	     
#---------------------------------------------------------------------------------------------------
# Compute Attributes Test	  

# Compute attribute test driver
set(${otb-module}ComputeAttributesTests
  	otbLSOBIAComputeAttributesTestDriver.cxx
    otbLSOBIAComputeAttributesTest.cxx
    )
    
add_executable(otbLSOBIAComputeAttributesTestDriver ${${otb-module}ComputeAttributesTests})
target_link_libraries(otbLSOBIAComputeAttributesTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIAComputeAttributesTestDriver)

otb_add_test(NAME otbLSOBIAComputeAttributesTest COMMAND otbLSOBIAComputeAttributesTestDriver
             otbLSOBIAComputeAttributesTest
             ${LSOBIA_SOURCE_DIR}/data/Vector_0_0.gml
             ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
             ${TEMP_DIR}
             ${TEMP_DIR}
             ${TEMP_DIR}/Vector_with_attributes_0_0.gml
)
