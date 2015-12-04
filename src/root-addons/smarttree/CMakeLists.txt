include_directories(.)

set(SRCS
    smarttree.h
  )

set(SRCS_EXTRA
    smarttree.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__tree.cc")
ROOT_GENERATE_DICTIONARY("${SRCS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})
set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++ -DR__ACCESS_IN_SYMBOL")
endif()

add_library(smarttree ${SRCS} ${SRCS_EXTRA})
target_link_libraries(smarttree base ${ROOT_LIBRARIES})
