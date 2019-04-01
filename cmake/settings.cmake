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

include(CheckCXXCompilerFlag)

# enable as many warnings as possible
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wnon-virtual-dtor -Werror -Wno-error=unused-variable")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers")  # fix compiling ROOT generated dictionaries with GCC >= 8
# cereal doesn't compile with GCC >= 8, temporarily fix this by ignoring the errors if the compiler supports this flag
# (remove once https://github.com/USCiLab/cereal/issues/497 is fixed)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
  CHECK_CXX_COMPILER_FLAG("-Wno-class-memaccess" COMPILER_SUPPORTS_CLASS_MEMACCESS)
  if(COMPILER_SUPPORTS_CLASS_MEMACCESS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-class-memaccess")
  endif()
endif()

# check which c++ standard is supported by the compiler
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
  # clang added warnings for unused lambda captures with the 5.0 release
  # the warnings are sometimes a bit too generous and warn about e.g. unused "this" even though its needed
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 5.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-lambda-capture")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=unused-but-set-variable")
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

# build shared libraries by default
SET(BUILD_SHARED_LIBS ON)

# disable the --as-needed globally, fixes linking problems on Ubuntu
set(DEFAULT_LINKER_FLAGS "-Wl,--no-as-needed")

# make sure the linker always finds ROOT's TDatabasePDG class
set(DEFAULT_LINKER_FLAGS "${DEFAULT_LINKER_FLAGS} -lEG")

# check if gold linker is available and use it
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -fuse-ld=gold -Wl,--version ERROR_QUIET OUTPUT_VARIABLE LD_VERSION)
if ("${LD_VERSION}" MATCHES "GNU gold")
    set(DEFAULT_LINKER_FLAGS "${DEFAULT_LINKER_FLAGS} -Wl,-fuse-ld=gold")
    message(STATUS "GNU google ld (gold) linker will be used.")
endif()
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${DEFAULT_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${DEFAULT_LINKER_FLAGS}")

# for file(GLOB_RECURSE..) don't follow symlinks
cmake_policy(SET CMP0009 NEW)

enable_testing()
# use some concurrency for tests
if(NOT CTEST_PARALLEL_JOBS)
  set(CTEST_PARALLEL_JOBS 2)
endif()
