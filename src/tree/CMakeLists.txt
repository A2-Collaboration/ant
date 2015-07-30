include_directories(.)

set(SRCS
  TDataRecord.h
  THeaderInfo.h
  TDetectorRead.h
  TSlowControl.h
  TUnpackerMessage.h
  TEvent.h
  TCluster.h
  TCandidate.h
  TTagger.h
  TTaggerHit.h
  TCalibrationData.h
  )

set(SRCS_EXTRA
  UnpackerReader.cc
  UnpackerWriter.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})
set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-DANT_TREE_ENABLE_COPY_CTOR -Wno-effc++")

add_library(tree ${SRCS} ${SRCS_EXTRA})
target_link_libraries(tree base ${ROOT_LIBRARIES})
