# every subdirectory has its own bin/lib path
# this should be changed to one "global" directory...
if(NOT DEFINED EXECUTABLE_OUTPUT_PATH)
	set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
endif()

if(NOT DEFINED LIBRARY_OUTPUT_PATH)
	set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
endif()

# we check for empty string here, since the variable
# is indeed defined to an empty string
if(NOT CMAKE_BUILD_TYPE)
  # this also reflects this default in the GUI
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# really no optimization in debug mode
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
else()
  message(FATAL_ERROR "Non-gnu compiler not supported at the moment")
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${BUILD_TYPE}})
