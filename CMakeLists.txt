cmake_minimum_required (VERSION 2.6)
project(ant)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

include(cmake/settings.cmake)

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# require a fairly recent ROOT version
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT "5.30" REQUIRED)
find_package(Pluto REQUIRED)
find_package(APLCON REQUIRED)


include_directories(src
  ${ROOT_INCLUDE_DIR}
  ${PLUTO_INCLUDE_DIR}
  ${APLCON_INCLUDE_DIR})

link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBS ${ROOT_LIBRARIES} ${APLCON_LIBRARY} ${PLUTO_LIBRARY})

add_subdirectory(src)
add_subdirectory(unpacker)