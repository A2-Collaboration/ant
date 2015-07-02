find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

set(SRCS
  main.cpp
  stl_helpers.h
  RawFileReader.cc
  RawFileReader.h
)

add_executable(unpacker ${SRCS})
target_link_libraries(unpacker ${LIBLZMA_LIBRARY})