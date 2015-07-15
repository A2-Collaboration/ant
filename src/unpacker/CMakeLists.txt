find_package(LibLZMA REQUIRED)
include_directories(${LIBLZMA_INCLUDE_DIRS})

include_directories(.)

set(SRCS_TREE
  tree/TDataRecord.h
  tree/THeaderInfo.h
  tree/TDetectorRead.h
  tree/TSlowControl.h
  tree/TUnpackerMessage.h
  tree/TEvent.h
  )

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

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS_TREE}" "${CMAKE_CURRENT_SOURCE_DIR}/tree/LinkDef.h" ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY} ${SRCS_TREE})

add_library(unpacker ${SRCS})
target_link_libraries(unpacker base expconfig ${LIBLZMA_LIBRARY} ${ROOT_LIBRARIES})

add_executable(unpacker_main main.cpp)
target_link_libraries(unpacker_main unpacker reconstruct)
