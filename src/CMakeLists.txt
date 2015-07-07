# this makes cross-includes easy,
# for example for like
# >>  #include "base/Logger.h"
# but don't forget to link against
# the corresponding library then

include_directories(.)

# the order is important if libraries from
# other subdirectories are needed
add_subdirectory(base)
add_subdirectory(expconfig)
add_subdirectory(unpacker)
add_subdirectory(analysis)
