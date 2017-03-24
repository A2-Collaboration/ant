# - Try to find APLCONpp
# Once done this will define
#  APLCONpp_FOUND - System has APLCONpp
#  APLCONpp_INCLUDE_DIRS - The APLCONpp include directories
#  APLCONpp_LIBRARIES - The libraries needed to use APLCONpp

set(APLCON_SEARCH_PATHS
  $ENV{APLCONSYS}
  ${CMAKE_SOURCE_DIR}/../APLCONpp)

find_path(APLCONpp_INCLUDE_DIR APLCON.hpp
  PATHS ${APLCON_SEARCH_PATHS}
  PATH_SUFFIXES src
  NO_DEFAULT_PATH
  )

find_library(APLCONpp_LIBRARY NAMES aplcon++
             PATHS ${APLCON_SEARCH_PATHS}
             PATH_SUFFIXES build
             NO_DEFAULT_PATH
            )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set APLCONpp_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(APLCONpp DEFAULT_MSG
                                  APLCONpp_LIBRARY APLCONpp_INCLUDE_DIR)

mark_as_advanced(APLCONpp_INCLUDE_DIR APLCONpp_LIBRARY)

set(APLCONpp_LIBRARIES ${APLCONpp_LIBRARY} )
set(APLCONpp_INCLUDE_DIRS ${APLCONpp_INCLUDE_DIR} )
