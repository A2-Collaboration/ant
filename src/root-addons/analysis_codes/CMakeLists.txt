include_directories(.)

set(DICT_HEADERS
    symmetric_pi0.h
    OmegaEtaG.h
    ExtractResolutions.h
    CBESum_Check.h
  )

set(SRCS
    symmetric_pi0.cc
    OmegaEtaG.cc
    TreeTools.cc
    ExtractResolutions.cc
    CBESum_Check.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__analysis_codes.cc")
ROOT_GENERATE_DICTIONARY("${DICT_HEADERS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})
set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++ -DR__ACCESS_IN_SYMBOL")
endif()

add_library(analysis_codes ${SRCS} ${DICT_HEADERS})
target_link_libraries(analysis_codes base analysis cbtaps_display ${ROOT_LIBRARIES})
