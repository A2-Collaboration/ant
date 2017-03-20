# every subdirectory has its own bin/lib path
# this should be changed to one "global" directory...
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

option(Ant_COVERAGE "Enable coverage build" OFF)

# we check for empty string here, since the variable
# is indeed defined to an empty string
if(Ant_COVERAGE)
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


# enable as many warnings as possible
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wnon-virtual-dtor -Werror")
# check which c++ standard is supported by the compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  message(STATUS "C++11 will be used: -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  message(STATUS "C++11 will be used: -std=c++0x")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# really no optimization in debug mode
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ftemplate-backtrace-limit=0")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # disable optimizations to fix clang infinite loop...
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
    message(STATUS "Disabling optimization for clang <3.8")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O1")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # for GCC >5.1 add -Wsuggest-override
else()
  message(FATAL_ERROR "Non-gnu compiler not supported at the moment")
endif()

option(Ant_MARCH "Enable auto-detection of CPU specific optimizations (march=native)" OFF)
if(Ant_MARCH)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${BUILD_TYPE}})

SET(BUILD_SHARED_LIBS ON)

# check if gold linker is available and use it
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -fuse-ld=gold -Wl,--version ERROR_QUIET OUTPUT_VARIABLE LD_VERSION)
if ("${LD_VERSION}" MATCHES "GNU gold")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-fuse-ld=gold")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-fuse-ld=gold")
    message(STATUS "GNU gold linker will be used.")
endif()

# for file(GLOB_RECURSE..) don't follow symlinks
cmake_policy(SET CMP0009 NEW)

enable_testing()
# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()
