find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

set(SRCS
  stl_helpers.h
  RawFileReader.cc
  RawFileReader.h
  Unpacker.cc
  Unpacker.h
  UnpackerAcquMk2.cc
  UnpackerAcquMk2.h
  TDataRecord.cc
  TDataRecord.h
)

include_directories(.)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker ${LIBLZMA_LIBRARY})

add_executable(unpacker_main main.cpp)
target_link_libraries(unpacker_main unpacker)