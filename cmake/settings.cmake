# every subdirectory has its own bin/lib path
# this should be changed to one "global" directory...
if(NOT DEFINED EXECUTABLE_OUTPUT_PATH)
        set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
endif()

if(NOT DEFINED LIBRARY_OUTPUT_PATH)
        set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
endif()

option(COVERAGE "Enable coverage build" OFF)

# we check for empty string here, since the variable
# is indeed defined to an empty string
if(COVERAGE)
  message(STATUS "Coverage build, enforce Debug build type")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
elseif(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()


# really no optimization in debug mode
if(CMAKE_COMPILER_IS_GNUCXX)
  # for GCC >5.1 add -Wsuggest-override
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
else()
  message(FATAL_ERROR "Non-gnu compiler not supported at the moment")
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${BUILD_TYPE}})

SET(BUILD_SHARED_LIBS ON)

# for file(GLOB_RECURSE..) don't follow symlinks
cmake_policy(SET CMP0009 NEW)

enable_testing()
# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()
