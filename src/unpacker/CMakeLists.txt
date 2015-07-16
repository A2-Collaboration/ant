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
  detail/UnpackerAcqu_detail.cc
  detail/UnpackerAcqu_detail.h
  detail/UnpackerAcqu_legacy.h
  detail/UnpackerAcqu_FileFormatMk1.cc
  detail/UnpackerAcqu_FileFormatMk1.h
  detail/UnpackerAcqu_FileFormatMk2.cc
  detail/UnpackerAcqu_FileFormatMk2.h
  detail/UnpackerAcqu_templates.h
)

add_library(unpacker ${SRCS})
target_link_libraries(unpacker base tree expconfig ${LIBLZMA_LIBRARY})

