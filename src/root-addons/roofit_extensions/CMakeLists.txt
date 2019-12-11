include_directories(.)

set(DICT_HEADERS
    RooGaussExp.h
    RooGaussDoubleSidedExp.h
  )

set(SRCS
    RooGaussExp.cxx
    RooGaussDoubleSidedExp.cxx
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__roofit_extensions.cc")
ROOT_GENERATE_DICTIONARY("${DICT_HEADERS}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
  ${ROOT_DICTIONARY})
list(APPEND SRCS ${ROOT_DICTIONARY})
set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++ -Wno-strict-aliasing")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${ROOT_DICTIONARY} PROPERTIES COMPILE_FLAGS
  "-Wno-effc++ -DR__ACCESS_IN_SYMBOL")
endif()


add_library(roofit_extensions ${SRCS} ${DICT_HEADERS})
target_link_libraries(roofit_extensions ${ROOT_LIBRARIES})
if (ROOT_VERSION VERSION_GREATER 6)
  add_custom_command(TARGET roofit_extensions POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/G__roofit_extensions_rdict.pcm ${PROJECT_BINARY_DIR}/lib
    )
  add_custom_command(TARGET roofit_extensions POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/libroofit_extensions.rootmap ${PROJECT_BINARY_DIR}/lib
    )
endif()
