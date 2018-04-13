include_directories(.)

set(DICT_HEADERS
    symmetric_pi0.h
    ExtractResolutions.h
    CBESum_Check.h
    TAPSEnergy_Check.h
    DetectorPlots.h
    BrowseHistograms.h
    PlotTimings.h
    hstack.h
    histtools.h
    ThetaCBToyMC.h
    KinfitExtract.h
    Math.h
    InterpolatedPulls.h
    DisplayClustering.h
    TwoPi0_MCSmearing_tools.h
    Omega.h
    pi0.h
    TimeDependentCalibration.h
    Fits.h
    pi0_true_calib.h
  )

set(SRCS
    symmetric_pi0.cc
    TreeTools.cc
    ExtractResolutions.cc
    CBESum_Check.cc
    TAPSEnergy_Check.cc
    DetectorPlots.cc
    BrowseHistograms.cc
    PlotTimings.cc
    hstack.cc
    histtools.cc
    ThetaCBToyMC.cc
    KinfitExtract.cc
    Math.cc
    InterpolatedPulls.cc
    DisplayClustering.cc
    TwoPi0_MCSmearing_tools.cc
    Omega.cc
    pi0.cc
    TimeDependentCalibration.cc
    Fits.cc
    hadd.cc
    pi0_true_calib.cc
  )

set(ROOT_DICTIONARY "${CMAKE_CURRENT_BINARY_DIR}/G__analysis_codes.cc")
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


add_library(analysis_codes ${SRCS} ${DICT_HEADERS})
target_link_libraries(analysis_codes base cbtaps_display reconstruct analysis_utils analysis_plot ${ROOT_LIBRARIES})
if (ROOT_VERSION VERSION_GREATER 6)
  add_custom_command(TARGET analysis_codes POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/G__analysis_codes_rdict.pcm ${PROJECT_BINARY_DIR}/lib
    )
endif()
