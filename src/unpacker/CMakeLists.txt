find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

include_directories(.)

set(SRCS
  RawFileReader.cc
  Unpacker.cc
  UnpackerA2Geant.cc
  UnpackerAcqu.cc
  detail/UnpackerAcqu_detail.cc
  detail/UnpackerAcqu_FileFormatMk1.cc
  detail/UnpackerAcqu_FileFormatMk2.cc
  detail/UnpackerAcqu_templates.h
  detail/UnpackerAcqu_legacy.h
)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker base tree expconfig ${LIBLZMA_LIBRARY})

