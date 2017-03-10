otb_module_test()


# ----------------------------------------------------------------------------------------------
# TODO : ADD: memory check tests (valgrind)
# TODO : ADD: otb_add_test_mpi (see OTBMPIConfig examples)
# TODO : ADD: use all the available drivers to add tests
# TODO : ADD: change app parameters to test both Baatz and MeanShift

otb_test_application(NAME apLSOBIA_BaatzOneProc
	APP LSSegmentation
	OPTIONS -io.im ${LSOBIA_SOURCE_DIR}/data/pleiades_1000_1000.tif
	-io.out ${TEMP}
	-io.temp ${TEMP}
	-algorithm baatz
	-algorithm.baatz.numitfirstpartial 5
	-algorithm.baatz.numitpartial 5
	-algorithm.baatz.stopping 40
	-algorithm.baatz.spectralweight 0.5
	-algorithm.baatz.geomweight 0.5
	-processing.memory 2000
	-processing.maxtilesizex 1000
	-processing.maxtilesizey 1000
	-processing.numproc 1)

# ----------------------------------------------------------------------------------------------
# Graph to image conversion test driver
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
	     ${TEMP}/image_from_graph.tif)

# ----------------------------------------------------------------------------------------------
# Large Scale Mean Shift test driver
set(${otb-module}LSMeanShiftTests
    otbLSOBIALSMeanShiftTestDriver.cxx
    otbObiaLSMeanShift.cxx
    )

add_executable(otbLSOBIALSMeanShiftTestDriver ${${otb-module}LSMeanShiftTests})
target_link_libraries(otbLSOBIALSMeanShiftTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIALSMeanShiftTestDriver)

# ----------------------------------------------------------------------------------------------
# Large Scale Small Regions Merging test driver
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
	     ${TEMP}/image_from_graph_without_gdal.tif)


# ----------------------------------------------------------------------------------------------
# Small Region Merging on a single file test driver
set(${otb-module}SmallRegionMergingOnFileTests
    otbLSOBIASmallRegionMergingOnFileTestDriver.cxx
    otbObiaSmallRegionMergingOnFile.cxx
    )

add_executable(otbLSOBIASmallRegionMergingOnFileTestDriver ${${otb-module}SmallRegionMergingOnFileTests})
target_link_libraries(otbLSOBIASmallRegionMergingOnFileTestDriver ${${otb-module}-Test_LIBRARIES})
otb_module_target_label(otbLSOBIASmallRegionMergingOnFileTestDriver)


# ----------------------------------------------------------------------------------------------
# Image Smoothing test driver
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



