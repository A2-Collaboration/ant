cmake_minimum_required (VERSION 2.8)
project(ant)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

include(cmake/settings.cmake)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

option(COVERALLS "Generate coveralls data" OFF)
enable_testing() # already here for correct coveralls functionality

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# require a fairly recent ROOT version
# this is needed by all subdirectories
find_package(ROOT "5.30" REQUIRED)
find_package(Pluto REQUIRED)
find_package(APLCON REQUIRED)

link_directories(${ROOT_LIBRARY_DIR})
include_directories(
  ${ROOT_INCLUDE_DIR}
  ${PLUTO_INCLUDE_DIR}
  ${APLCON_INCLUDE_DIR}
)

# now the relevant subdirectories
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(progs)

include(cmake/doxygen.cmake)

# create a rootlogon.C
configure_file(cmake/rootlogon.C.in rootlogon.C)

