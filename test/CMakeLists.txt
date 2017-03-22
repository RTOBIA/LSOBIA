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
	-io.out ${TEMP_DIR}
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
	-io.out ${TEMP_DIR}
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


# ----------------------------------------------------------------------------------------------
# Segmentation Application Tests
# - Mean Shift Tests
otb_test_application(NAME apLSOBIA_MeanShiftOneProc
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out ${TEMP_DIR}
	-io.temp ${TEMP_DIR}
	-algorithm meanshift
	-algorithm.meanshift.maxiter 100
	-algorithm.meanshift.spatialr 15
	-algorithm.meanshift.spectralr 16
	-algorithm.meanshift.threshold 0.01
	-algorithm.meanshift.ranger 0
	-algorithm.meanshift.modesearch "off"
	-processing.memory 2000
<<<<<<< HEAD
	-processing.maxtilesizex 200
	-processing.maxtilesizey 200
	-processing.numproc 1)
=======
	-processing.maxtilesizex 1000
	-processing.maxtilesizey 1000
	)
>>>>>>> 9c412f368f0e60a0c5cdf0ea98db137cb44c294f

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
# Image Smoothing test driver
# TODO : implement driver and add tests
set(${otb-module}SmoothImageTests
    otbLSOBIASmoothImageTestDriver.cxx
    otbObiaSmoothImage.cxx
    )

add_executable(otbLSOBIASmoothImageTestDriver ${${otb-module}SmoothImageTests})
target_link_libraries(otbLSOBIASmoothImageTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIASmoothImageTestDriver)


# ----------------------------------------------------------------------------------------------
# Graph Splitting test driver
set(${otb-module}SplitGraphTests
    otbLSOBIASplitGraphTestDriver.cxx
    otbObiaSplitGraph.cxx
    )


add_executable(otbLSOBIASplitGraphTestDriver ${${otb-module}SplitGraphTests})
target_link_libraries(otbLSOBIASplitGraphTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIASplitGraphTestDriver)

<<<<<<< HEAD

# ----------------------------------------------------------------------------------------------
# Small Regions Merging Test

otb_test_application(NAME apLSOBIA_SRMOneProc
	APP LSSmallRegionsMerging
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out ${TEMP_DIR}
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
	-processing.numproc 1)




=======
otb_add_test(NAME otbLSOBIASplitGraph COMMAND otbLSOBIASplitGraphTestDriver
             otbObiaSplitGraph
             ${LSOBIA_SOURCE_DIR}/data/Graph_0_0.dat
             ${TEMP_DIR}
	     250
             250
             )
>>>>>>> 9c412f368f0e60a0c5cdf0ea98db137cb44c294f
