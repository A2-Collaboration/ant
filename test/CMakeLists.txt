# the framework for our tests
include_directories(.)

add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND} -V)
add_library(catch main.cpp catch.hpp)
enable_testing()

# needs file and name for test
set(TESTNAME unpacker_RawFileReader)
add_executable(${TESTNAME} unpacker/RawFileReader.cc)
target_link_libraries(${TESTNAME} catch)
add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
add_dependencies(build_and_test ${TESTNAME})