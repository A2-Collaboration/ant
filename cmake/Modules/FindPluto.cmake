# - Try to find GSI Hades Pluto
# Once done this will define
#  PLUTO_FOUND - System has Pluto
#  PLUTO_INCLUDE_DIRS - The Pluto include directories
#  PLUTO_LIBRARIES - The libraries needed to use Pluto

Set(PLUTO_SEARCHPATH
        $ENV{PLUTOSYS}
        $ENV{HOME}/src/pluto
        $ENV{HOME}/opt/pluto
        $ENV{HOME}/pluto
        /opt/pluto/
)

find_library(PLUTO_LIBRARY NAMES Pluto PATHS ${PLUTO_SEARCHPATH})
find_path(PLUTO_INCLUDE_DIR NAMES PPlutoBulkDecay.h PATHS ${PLUTO_SEARCHPATH} PATH_SUFFIXES src include)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PLUTO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Pluto DEFAULT_MSG
                                  PLUTO_LIBRARY PLUTO_INCLUDE_DIR)

mark_as_advanced(PLUTO_INCLUDE_DIR PLUTO_LIBRARY)

set(PLUTO_LIBRARIES ${PLUTO_LIBRARY} )
set(PLUTO_INCLUDE_DIRS ${PLUTO_INCLUDE_DIR} )
