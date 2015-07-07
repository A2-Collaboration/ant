find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

include_directories(.)

set(SRCS
  RawFileReader.cc
  RawFileReader.h
  Unpacker.cc
  Unpacker.h
  UnpackerAcqu.cc
  UnpackerAcqu.h
  detail/UnpackerAcqu.cc
  detail/UnpackerAcqu.h
  detail/UnpackerAcqu_legacy.h
  tree/TDataRecord.h
  tree/THeaderInfo.h
  tree/TDetectorRead.h
  tree/TSlowControl.h
  tree/TUnpackerMessage.h
)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker base expconfig ${LIBLZMA_LIBRARY})

add_executable(unpacker_main main.cpp)
target_link_libraries(unpacker_main unpacker)
