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

set(ANT_BASE
  src/base/types.cc
  src/base/types.h
  src/ParticleType.h
  src/ParticleType.cc
  src/data/Track.h
  src/data/Track.cc
  src/data/Particle.h
  src/data/Particle.cc
  src/data/Event.h
  src/data/Event.cc
  src/Physics.h
  src/Physics.cc
  src/base/printable.h
  src/base/printable.cc
  src/utils/combinatorics.h
  src/plot/Histogram.h
  src/plot/Histogram.cc
  src/base/interval.h
  src/data/TaggerHit.h
  src/data/TaggerHit.cc
  src/plot/root_draw.h
  src/plot/root_draw.cc
  src/Detector.h
  src/Detector.cc
  src/utils/matcher.h
  src/data/TriggerInfo.h
  src/data/TriggerInfo.cc
  src/A2GeoAcceptance.h
  src/A2GeoAcceptance.cc
  src/plot/SmartHist.h
  src/plot/SmartHist.cc
  src/plot/HistogramFactories.h
  src/plot/HistogramFactories.cc
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
  src/input/DataReader.h
  src/input//DataReader.cc
  src/OutputManager.h
  src/OutputManager.cc
  )

set(ANT_PHYSICS
  src/analysis/TestPhysics.h
  src/analysis/TestPhysics.cc
  src/analysis/DeltaPlusPhysics.h
  src/analysis/DeltaPlusPhysics.cc
  src/analysis/omega.h
  src/analysis/omega.cc
  src/analysis/omega_bottomup.h
  src/analysis/omega_bottomup.cc
  src/analysis/GeoAcceptance.h
  src/analysis/GeoAcceptance.cc
  src/analysis/TestAPLCON.h
  src/analysis/TestAPLCON.cc
  src/analysis/DataOverview.h
  src/analysis/DataOverview.cc
  src/analysis/mctrue_acceptance.h
  src/analysis/mctrue_acceptance.cc
  )

add_library(antlib SHARED ${ANT_BASE} ${ANT_PHYSICS})
list(APPEND LIBS antlib)

aux_source_directory("src/progs" ANT_PROGS)
foreach(PROG_SRC ${ANT_PROGS})

  get_filename_component(name ${PROG_SRC} NAME_WE)
  message(STATUS "Adding ant executable: ${name}")
  add_executable(${name}
    ${PROG_SRC}
    )

  target_link_libraries(${name} ${LIBS})
endforeach()
