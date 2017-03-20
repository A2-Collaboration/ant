# the catch++ library is the workhorse for our tests
include_directories(.)
add_custom_target(build_and_test
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -j${CTEST_PARALLEL_JOBS}
  COMMENT "Running build_and_test with ${CTEST_PARALLEL_JOBS} workers"
  )
add_library(catch EXCLUDE_FROM_ALL catch.cc catch_config.h)
target_link_libraries(catch base third_party)

add_library(expconfig_helpers EXCLUDE_FROM_ALL expconfig_helpers.cc)
target_link_libraries(expconfig_helpers expconfig)

# some tests need binary blobs
# use a configure file to
set(TEST_BLOBS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/_blobs")
configure_file(catch_config.h.in catch_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# two helpers for defining tests
# we don't simply GLOB over all files,
# because macros are more flexible

macro(add_test_subdirectory subdir)
  set(TEST_SUBDIR ${subdir})
  include_directories(${CMAKE_SOURCE_DIR}/src/${subdir})
  add_subdirectory(${subdir})
endmacro()

macro(add_ant_test name)
  set(TESTDIR "${CMAKE_BINARY_DIR}/bin_test")
  set(TESTNAME "${TEST_SUBDIR}_${name}")
  set(TESTTARGET "test_${TESTNAME}")
  set(TESTFILE "Test${name}.cc")
  add_executable(${TESTTARGET} EXCLUDE_FROM_ALL ${TESTFILE})
  target_link_libraries(${TESTTARGET} catch ${TEST_SUBDIR} ${ARGN})
  set_target_properties(${TESTTARGET}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${TESTDIR}
    OUTPUT_NAME ${TESTNAME}
  )
  # search for expconfig, then automatically add expconfig_helpers
  set(_list ${TEST_SUBDIR} ${ARGN})
  list(FIND _list "expconfig" _index)
  if(${_index} GREATER -1)
      target_link_libraries(${TESTTARGET} expconfig_helpers)
  endif()
  add_test(NAME ${TESTTARGET} COMMAND "${TESTDIR}/${TESTNAME}")
  set_tests_properties(${TESTTARGET} PROPERTIES
    TIMEOUT 300
    WORKING_DIRECTORY ${TESTDIR})
  add_dependencies(build_and_test ${TESTTARGET})
endmacro()

macro(add_python_test_directory name)
  set(TESTNAME "test_python_${name}")
  find_program(PYCOV coverage)
  if(NOT PYCOV AND NOT PYTHONINTERP_FOUND)
    message(STATUS "Python interpreter not found. Skip Python tests.")
  elseif(NOT PYCOV)
    add_test(NAME ${TESTNAME} COMMAND ${PROJECT_SOURCE_DIR}/test/test_python_helper.sh ${name} ${PROJECT_SOURCE_DIR})
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 60)
  else()
    add_test(NAME ${TESTNAME} COMMAND ${PROJECT_SOURCE_DIR}/test/test_python_helper.sh ${name} ${PROJECT_SOURCE_DIR} 1)
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 60)
  endif()
endmacro()

include_directories(${CMAKE_SOURCE_DIR}/src)

add_test_subdirectory(base)
add_test_subdirectory(unpacker)
add_test_subdirectory(tree)
add_test_subdirectory(expconfig)
add_test_subdirectory(reconstruct)
add_test_subdirectory(analysis)
add_test_subdirectory(calibration)
add_python_test_directory(extra)

