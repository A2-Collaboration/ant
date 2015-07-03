# the framework for our tests
include(GetCatch)
include_directories(${CATCH_INCLUDE_DIR})

add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND} -V)
add_library(catch_main main.cpp)
add_dependencies(catch_main catch)
enable_testing()

# needs file and name for test
set(TESTNAME unpacker_RawFileReader)
add_executable(${TESTNAME} unpacker/RawFileReader.cc)
target_link_libraries(${TESTNAME} catch_main)
add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
add_dependencies(build_and_test ${TESTNAME})