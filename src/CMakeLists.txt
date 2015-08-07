# this makes cross-includes easy,
# for example for like
# >>  #include "base/Logger.h"
# but don't forget to link against
# the corresponding library then
include_directories(.)

# setup coverage
if(COVERALLS)
  message(STATUS "Code coverage testing enabling")
  include(Coveralls)
  coveralls_turn_on_coverage()   
endif()

macro(add_ant_subdirectory name)
  message(STATUS "Configuring ${name}...")
  add_subdirectory(${name})
  if(COVERALLS)
    # in base we exclude the subdirectories
    if(${name} STREQUAL "base")
      file(GLOB SRCS "${name}/*.cc" "${name}/*.h")
      list(APPEND COVERAGE_SRCS ${SRCS})
    else()
      file(GLOB_RECURSE SRCS "${name}/*.cc" "${name}/*.h")
      list(APPEND COVERAGE_SRCS ${SRCS})
    endif()
  endif()
endmacro()

# the order is important if libraries from
# other subdirectories are needed
add_ant_subdirectory(base)
add_ant_subdirectory(expconfig)
add_ant_subdirectory(tree)
add_ant_subdirectory(rint)
add_ant_subdirectory(unpacker)
add_ant_subdirectory(reconstruct)
add_ant_subdirectory(analysis)
add_ant_subdirectory(calibration)


if(COVERALLS)
  # Create the coveralls target.
  coveralls_setup(
    "${COVERAGE_SRCS}" # The source files.
    ON                 # If we should upload.
    )   
endif()
