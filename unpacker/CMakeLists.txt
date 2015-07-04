find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

set(SRCS
  stl_helpers.h
  RawFileReader.cc
  RawFileReader.h
  Unpacker.cc
  Unpacker.h
  UnpackerAcqu.cc
  UnpackerAcqu.h
  UnpackerAcqu_detail.h # move in subfolder later
  TDataRecord.cc
  TDataRecord.h
)

include_directories(.)

add_library(logger
  Logger.cc
  Logger.h
  easylogging++.h
)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker logger ${LIBLZMA_LIBRARY})

add_executable(unpacker_main main.cpp)
target_link_libraries(unpacker_main unpacker)
