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
  UnpackerAcqu_detail.cc
  UnpackerAcqu_detail.h
  UnpackerAcqu_legacy.h
  TDataRecord.cc
  TDataRecord.h
)

include_directories(.)

# hack in some useful ant base stuff for now
include_directories(../src/base)

add_library(logger
  ../src/base/printable.cc
  Logger.cc
  Logger.h
  easylogging++.h
)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker logger ${LIBLZMA_LIBRARY})

add_executable(unpacker_main main.cpp)
target_link_libraries(unpacker_main unpacker)
