cmake_minimum_required(VERSION 2.8.9)
project(LSOBIA)

list(APPEND CMAKE_MODULE_PATH ${LSOBIA_SOURCE_DIR})

set(LSOBIA_LIBRARIES LSOBIA)

if(NOT OTB_SOURCE_DIR)
    find_package(OTB REQUIRED)
    include(OTBModuleExternal)
else()
    otb_module_impl()
endif()
