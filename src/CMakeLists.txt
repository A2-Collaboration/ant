# this makes cross-includes easy,
# for example for like
# >>  #include "base/Logger.h"
# but don't forget to link against
# the corresponding library then
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# for generated files by configure_file()
include_directories(${CMAKE_CURRENT_BINARY_DIR})


macro(add_ant_subdirectory name)
  message(STATUS "Configuring ${name}...")
  add_subdirectory(${name})
endmacro()

# the order is important if libraries from
# other subdirectories are needed
add_ant_subdirectory(base)
add_ant_subdirectory(expconfig)
add_ant_subdirectory(tree)
add_ant_subdirectory(unpacker)
add_ant_subdirectory(reconstruct)
add_ant_subdirectory(analysis)
add_ant_subdirectory(calibration)
add_ant_subdirectory(simulation/mc)
add_ant_subdirectory(smarttree)
