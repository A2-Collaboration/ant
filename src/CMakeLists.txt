# this makes cross-includes easy,
# for example for like
# >>  #include "base/Logger.h"
# but don't forget to link against
# the corresponding library then
include_directories(.)

# require a fairly recent ROOT version
# this is needed by all subdirectories
find_package(ROOT "5.30" REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})


# the order is important if libraries from
# other subdirectories are needed
add_subdirectory(base)
add_subdirectory(expconfig)
add_subdirectory(unpacker)
add_subdirectory(analysis)
