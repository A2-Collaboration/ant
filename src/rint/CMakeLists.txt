include_directories(.)

set(SRCS
    TCalibrationEditor.h
  )

set(SRCS_EXTRA
    TCalibrationEditor.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__rint.cc")
ROOT_GENERATE_DICTIONARY("${SRCS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})
set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-DANT_TREE_ENABLE_COPY_CTOR -Wno-effc++")

add_library(rint ${SRCS} ${SRCS_EXTRA})
target_link_libraries(rint base calibration ${ROOT_LIBRARIES})
