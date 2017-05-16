include_directories(.)

set(SRCS_DICT
  TID.h
  TEvent.h
  TCalibrationData.h
  TAntHeader.h
  ${ROOT_INCLUDE_DIR}/TLorentzVector.h
  )

set(SRCS
  MemoryPool.h
  stream_TBuffer.h
  TDetectorReadHit.h
  TSlowControl.h
  TUnpackerMessage.h
  TTarget.h
  TTrigger.h
  TTaggerHit.h
  TCluster.h
  TCandidate.h
  TParticle.cc
  TSimpleParticle.cc
  TEventData.cc
  TEvent.cc
  TAntHeader.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS_DICT}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})

set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++ -DR__ACCESS_IN_SYMBOL")
endif()


add_library(tree ${SRCS} ${SRCS_DICT})
target_link_libraries(tree base ${ROOT_LIBRARIES})
