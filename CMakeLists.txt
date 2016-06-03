cmake_minimum_required(VERSION 2.8)
project( CompressedPolyDataServer )
get_filename_component(PARENT_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
set(VTK_DIR "${PARENT_DIR}/Slicer-unixBuild/VTKv6-build")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

get_filename_component(PARENT_DIR ${CMAKE_SOURCE_DIR} DIRECTORY)
message("${PARENT_DIR}/OpenIGTLink-XcodeBuild/Testing/OpenH264")
set(OpenH264_INCLUDE "${PARENT_DIR}/OpenIGTLink-XcodeBuild/Testing/OpenH264" CACHE STRING "OpenH264 source file")
set(OpenH264_LIB "${PARENT_DIR}/OpenIGTLink-XcodeBuild/Testing/OpenH264/Debug" CACHE STRING "OpenH264 library")

include_directories("${OpenH264_INCLUDE}/codec")
include_directories("${OpenH264_INCLUDE}/res")
include_directories("${OpenH264_INCLUDE}/test")
LINK_DIRECTORIES("${OpenH264_LIB}")

set(OpenIGTLink_INCLUDE "${PARENT_DIR}/OpenIGTLink/Source" CACHE STRING "OpenIGTLink source file")
set(OpenIGTLink_LIB "${PARENT_DIR}/OpenIGTLink-XcodeBuild" CACHE STRING "OpenIGTLink library")

include_directories("${OpenIGTLink_INCLUDE}")
include_directories("${OpenIGTLink_INCLUDE}/igtlutil")
include_directories("${OpenIGTLink_LIB}")
LINK_DIRECTORIES("${OpenIGTLink_LIB}/bin")

set(pcl_common_INCLUDE "/usr/local/include/pcl-1.6/" CACHE STRING "plc_common source file")
set(plc_common_LIB "/usr/local/lib" CACHE STRING "plc_common library")
include_directories("${pcl_common_INCLUDE}")
LINK_DIRECTORIES("${plc_common_LIB}")

set(EIGEN_INCLUDE "/opt/local/include/eigen3/" CACHE STRING "EIGEN source file")
set(EIGEN_LIB "/opt/local/lib" CACHE STRING "EIGEN library")
include_directories("${EIGEN_INCLUDE}")
LINK_DIRECTORIES("${EIGEN_LIB}")

set(boost_INCLUDE "/opt/local/include/" CACHE STRING "boost source file")
set(boost_LIB "/opt/local/lib" CACHE STRING "boost library")
include_directories("${boost_INCLUDE}")
LINK_DIRECTORIES("${boost_LIB}")

add_library(plyReader vtkPLYReaderLocal.cpp)
add_library(range_image_local range_image_local.cpp)
add_executable( CompressedPolyDataServer CompressedPolyDataServer )
target_link_libraries( CompressedPolyDataServer OpenIGTLink openH264 plyReader ${plc_common_LIB} range_image_local)
