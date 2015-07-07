# the framework for our tests
include_directories(.)
add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND} -V)
add_library(catch main.cpp catch.hpp)
enable_testing()

# two helpers for defining tests
# we don't simply GLOB over all files,
# because macros are more flexible

macro(add_test_subdirectory subdir)
  set(TEST_SUBDIR ${subdir})
  include_directories(${CMAKE_SOURCE_DIR}/src/${subdir})
  add_subdirectory(${subdir})
endmacro()

macro(add_ant_test name)
  set(TESTNAME "test_${TEST_SUBDIR}_${name}")
  set(TESTFILE "Test${name}.cc")
  add_executable(${TESTNAME} ${TESTFILE})
  target_link_libraries(${TESTNAME} catch ${TEST_SUBDIR})
  add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 5)
  add_dependencies(build_and_test ${TESTNAME})
endmacro()

add_test_subdirectory(unpacker)

