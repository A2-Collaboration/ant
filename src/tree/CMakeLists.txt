include_directories(.)

set(SRCS
  TDataRecord.h
  THeaderInfo.h
  TDetectorRead.h
  TSlowControl.h
  TUnpackerMessage.h
  TEvent.h
  TCluster.h
  TTrack.h
  TTaggerHit.h
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})

add_library(tree ${SRCS})
target_link_libraries(tree base ${ROOT_LIBRARIES})
