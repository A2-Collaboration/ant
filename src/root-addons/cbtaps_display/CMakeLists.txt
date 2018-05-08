include_directories(.)


set(SRCS
    TH2Crystals.h
    TH2CB.h
    TH2TAPS.h
)

set(SRCS_EXTRA
    TH2TAPS.cc
    TH2DrawTool.cc
    TH2CB.cc
    TH2Crystals.cc
)


set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__cbtaps_display.cc")
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


add_library(cbtaps_display ${SRCS} ${SRCS_EXTRA})
target_link_libraries(cbtaps_display base ${ROOT_LIBRARIES})
if (ROOT_VERSION VERSION_GREATER 6)
  add_custom_command(TARGET cbtaps_display POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/G__cbtaps_display_rdict.pcm ${PROJECT_BINARY_DIR}/lib
    )
  add_custom_command(TARGET cbtaps_display POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/libcbtaps_display.rootmap ${PROJECT_BINARY_DIR}/lib
    )
endif()
