cmake_minimum_required (VERSION 3.0)
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

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Linker Flags:    " ${CMAKE_SHARED_LINKER_FLAGS})
message(STATUS "*** Install libs to: " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "*** Install bin to:  " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# require a fairly recent ROOT version
# this is needed by all subdirectories
find_package(ROOT "5.30" REQUIRED)
find_package(Pluto REQUIRED)
find_package(APLCONpp REQUIRED)
find_package(GSL REQUIRED)

include(cmake/check_dependencies.cmake)

link_directories(${ROOT_LIBRARY_DIR})
# including them as SYSTEM prevents
# many false-positive warnings
include_directories(SYSTEM
  ${ROOT_INCLUDE_DIR}
  ${PLUTO_INCLUDE_DIR}
  ${APLCONpp_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
)

# now the relevant subdirectories
add_subdirectory(third-party)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(progs)

include(cmake/doxygen.cmake)

# create a rootlogon.C
configure_file(cmake/rootlogon.C.in rootlogon.C)

# link the AntSubmit file
add_custom_target(link_AntSubmit ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  "${CMAKE_SOURCE_DIR}/extra/AntSubmit"
  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/AntSubmit"
  )

# link the AntSimSubmit file
add_custom_target(link_AntSimSubmit ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  "${CMAKE_SOURCE_DIR}/extra/SimulationBlaster/Ant_simulation_blaster.py"
  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/AntSimSubmit"
  )

# link the AntMapReduce file
add_custom_target(link_AntMapReduce ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  "${CMAKE_SOURCE_DIR}/extra/AntMapReduce"
  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/AntMapReduce"
  )
