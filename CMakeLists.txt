cmake_minimum_required (VERSION 2.6)
project(ant)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

include(cmake/settings.cmake)

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# require a fairly recent ROOT version
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT "5.30" REQUIRED)
find_package(Pluto REQUIRED)
find_package(APLCON REQUIRED)


include_directories(src
  ${ROOT_INCLUDE_DIR}
  ${PLUTO_INCLUDE_DIR}
  ${APLCON_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBS ${ROOT_LIBRARIES} ${APLCON_LIBRARY} ${PLUTO_LIBRARY})

set(ANALYSIS_BASE
  src/base/types.cc
  src/base/types.h
  src/base/interval.h
  src/base/printable.h
  src/base/printable.cc
  src/ParticleType.h
  src/ParticleType.cc
  src/physics/Physics.h
  src/physics/Physics.cc
  src/utils/combinatorics.h
  src/utils/matcher.h
  src/Detector.h
  src/Detector.cc
  src/A2GeoAcceptance.h
  src/A2GeoAcceptance.cc
  src/input/DataReader.h
  src/input/DataReader.cc
  src/OutputManager.h
  src/OutputManager.cc
  )

set(ANALYSIS_PLOT
  src/plot/Histogram.h
  src/plot/Histogram.cc
  src/plot/root_draw.h
  src/plot/root_draw.cc
  src/plot/SmartHist.h
  src/plot/SmartHist.cc
  src/plot/HistogramFactories.h
  src/plot/HistogramFactories.cc
)

set(ANALYSIS_DATASTRUCTURES
  src/data/Track.h
  src/data/Track.cc
  src/data/Particle.h
  src/data/Particle.cc
  src/data/Event.h
  src/data/Event.cc
  src/data/TaggerHit.h
  src/data/TaggerHit.cc
  src/data/TriggerInfo.h
  src/data/TriggerInfo.cc
)

set(ANALYSIS_INPUT_GOAT
  src/input/goat/detail/InputModule.h
  src/input/goat/detail/InputModule.cc
  src/input/goat/detail/FileManager.h
  src/input/goat/detail/FileManager.cc
  src/input/goat/detail/PlutoInput.h
  src/input/goat/detail/PlutoInput.cc
  src/input/goat/detail/TreeManager.h
  src/input/goat/detail//TreeManager.cc
  src/input/goat/detail/TrackInput.h
  src/input/goat/detail/TrackInput.cc
  src/input/goat/detail/ParticleInput.h
  src/input/goat/detail/ParticleInput.cc
  src/input/goat/detail/TriggerInput.h
  src/input/goat/detail/TriggerInput.cc
  src/input/goat/detail/DetectorHitInput.h
  src/input/goat/detail/DetectorHitInput.cc
  src/input/goat/detail/GeantInput.h
  src/input/goat/detail/GeantInput.cc
  src/input/goat/detail/TaggerInput.h
  src/input/goat/detail/TaggerInput.cc
  src/input/goat/GoatReader.h
  src/input/goat/GoatReader.cc
)

set(ANALYSIS_PHYSICS
  src/physics/test/TestPhysics.h
  src/physics/test/TestPhysics.cc
  src/physics/DeltaPlus/DeltaPlusPhysics.h
  src/physics/DeltaPlus/DeltaPlusPhysics.cc
  src/physics/omega/omega.h
  src/physics/omega/omega.cc
  src/physics/omega/omega_bottomup.h
  src/physics/omega/omega_bottomup.cc
  src/physics/omega/GeoAcceptance.h
  src/physics/omega/GeoAcceptance.cc
  src/physics/test/TestAPLCON.h
  src/physics/test/TestAPLCON.cc
  src/physics/common/DataOverview.h
  src/physics/common/DataOverview.cc
  src/physics/omega/mctrue_acceptance.h
  src/physics/omega/mctrue_acceptance.cc
  )

set(ANALYSIS_ALL
    ${ANALYSIS_BASE}
    ${ANALYSIS_PLOT}
    ${ANALYSIS_DATASTRUCTURES}
    ${ANALYSIS_INPUT_GOAT}
    ${ANALYSIS_PHYSICS}
)

add_library(analysis SHARED ${ANALYSIS_ALL})
list(APPEND LIBS analysis)

aux_source_directory("src/progs" ANALYSIS_PROGS)
foreach(PROG_SRC ${ANALYSIS_PROGS})

  get_filename_component(name ${PROG_SRC} NAME_WE)
  message(STATUS "Adding ant executable: ${name}")
  add_executable(${name}
    ${PROG_SRC}
    )

  target_link_libraries(${name} ${LIBS})
endforeach()
