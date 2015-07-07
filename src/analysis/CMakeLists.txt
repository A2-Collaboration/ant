# analysis needs some more libraries
find_package(Pluto REQUIRED)
find_package(APLCON REQUIRED)

include_directories(
  .
  ${PLUTO_INCLUDE_DIR}
  ${APLCON_INCLUDE_DIR}
  )


set(ANALYSIS_BASE
  ParticleType.h
  ParticleType.cc
  physics/Physics.h
  physics/Physics.cc
  utils/combinatorics.h
  utils/matcher.h
  Detector.h
  Detector.cc
  A2GeoAcceptance.h
  A2GeoAcceptance.cc
  input/DataReader.h
  input/DataReader.cc
  OutputManager.h
  OutputManager.cc
  )

set(ANALYSIS_PLOT
  plot/Histogram.h
  plot/Histogram.cc
  plot/root_draw.h
  plot/root_draw.cc
  plot/SmartHist.h
  plot/SmartHist.cc
  plot/HistogramFactories.h
  plot/HistogramFactories.cc
)

set(ANALYSIS_DATASTRUCTURES
  data/Track.h
  data/Track.cc
  data/Particle.h
  data/Particle.cc
  data/Event.h
  data/Event.cc
  data/TaggerHit.h
  data/TaggerHit.cc
  data/TriggerInfo.h
  data/TriggerInfo.cc
)

set(ANALYSIS_INPUT_GOAT
  input/goat/detail/InputModule.h
  input/goat/detail/InputModule.cc
  input/goat/detail/FileManager.h
  input/goat/detail/FileManager.cc
  input/goat/detail/PlutoInput.h
  input/goat/detail/PlutoInput.cc
  input/goat/detail/TreeManager.h
  input/goat/detail//TreeManager.cc
  input/goat/detail/TrackInput.h
  input/goat/detail/TrackInput.cc
  input/goat/detail/ParticleInput.h
  input/goat/detail/ParticleInput.cc
  input/goat/detail/TriggerInput.h
  input/goat/detail/TriggerInput.cc
  input/goat/detail/DetectorHitInput.h
  input/goat/detail/DetectorHitInput.cc
  input/goat/detail/GeantInput.h
  input/goat/detail/GeantInput.cc
  input/goat/detail/TaggerInput.h
  input/goat/detail/TaggerInput.cc
  input/goat/GoatReader.h
  input/goat/GoatReader.cc
)

set(ANALYSIS_PHYSICS
  physics/test/TestPhysics.h
  physics/test/TestPhysics.cc
  physics/DeltaPlus/DeltaPlusPhysics.h
  physics/DeltaPlus/DeltaPlusPhysics.cc
  physics/omega/omega.h
  physics/omega/omega.cc
  physics/omega/omega_bottomup.h
  physics/omega/omega_bottomup.cc
  physics/omega/GeoAcceptance.h
  physics/omega/GeoAcceptance.cc
  physics/test/TestAPLCON.h
  physics/test/TestAPLCON.cc
  physics/common/DataOverview.h
  physics/common/DataOverview.cc
  physics/omega/mctrue_acceptance.h
  physics/omega/mctrue_acceptance.cc
  )

set(ANALYSIS_ALL
    ${ANALYSIS_BASE}
    ${ANALYSIS_PLOT}
    ${ANALYSIS_DATASTRUCTURES}
    ${ANALYSIS_INPUT_GOAT}
    ${ANALYSIS_PHYSICS}
)

add_library(analysis ${ANALYSIS_ALL})
target_link_libraries(analysis
  ${ROOT_LIBRARIES}
  ${PLUTO_LIBRARY}
  ${APLCON_LIBRARY}  
  base
  )

aux_source_directory("progs" ANALYSIS_PROGS)
foreach(PROG_SRC ${ANALYSIS_PROGS})

  get_filename_component(name ${PROG_SRC} NAME_WE)
  message(STATUS "Adding ant executable: ${name}")
  add_executable(${name}
    ${PROG_SRC}
    )

  target_link_libraries(${name} analysis)
endforeach()
