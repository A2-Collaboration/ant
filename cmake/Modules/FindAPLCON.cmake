# find the APLCON++ constrained fitter wrapper
# simplistic, should mark variables as cached or advanced

message(STATUS "Looking for APLCON...")
set(APLCON_FOUND FALSE)

set(APLCON_SEARCH_PATHS
  ${CMAKE_SOURCE_DIR}/../APLCON)

find_path(APLCON_INCLUDE_DIR APLCON.hpp
  PATHS ${APLCON_SEARCH_PATHS}/src
  NO_DEFAULT_PATH
  )
if(NOT APLCON_INCLUDE_DIR)
  Message(STATUS "Looking for APLCON... - APLCON.hpp not found")
  return()
endif()

find_library(APLCON_LIBRARY NAMES aplcon++
  PATHS ${APLCON_SEARCH_PATHS}
  PATH_SUFFIXES build
  NO_DEFAULT_PATH
  )
if(NOT APLCON_LIBRARY)
  Message(STATUS "Looking for APLCON... - aplcon++ library not found.")
endif()


if(APLCON_INCLUDE_DIR AND APLCON_LIBRARY)
  set(APLCON_FOUND TRUE)
  Message(STATUS "Looking for APLCON... - Found ${APLCON_LIBRARIES}")
endif()

if(NOT APLCON_FOUND AND APLCON_FIND_REQUIRED)
  message(FATAL_ERROR "APLCON is required, exit.")
endif()