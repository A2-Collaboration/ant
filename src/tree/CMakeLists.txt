include_directories(.)

set(SRCS
  TID.h
  TDataRecord.h
  TDetectorRead.h
  TSlowControl.h
  TUnpackerMessage.h
  TEvent.h
  TTagger.h
  TTaggerHit.h
  TCalibrationData.h
  TAntHeader.h
  ${ROOT_INCLUDE_DIR}/TLorentzVector.h
  )

set(SRCS_EXTRA
  MemoryPool.h
  UnpackerReader.cc
  UnpackerWriter.cc
  TCluster.h
  TCandidate.h
  TParticle.cc
  TEventData.h
  TEvent.cc
  TAntHeader.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})

set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-DANT_TREE_ENABLE_COPY_CTOR -Wno-effc++")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-DANT_TREE_ENABLE_COPY_CTOR -Wno-effc++ -DR__ACCESS_IN_SYMBOL")
endif()

add_library(tree ${SRCS} ${SRCS_EXTRA})
target_link_libraries(tree base ${ROOT_LIBRARIES})
