include_directories(
  .
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
  data/Candidate.h
  data/Candidate.cc
  data/Particle.h
  data/Particle.cc
  data/Event.h
  data/Event.cc
  data/TaggerHit.h
  data/TaggerHit.cc
  data/TriggerInfo.h
  data/TriggerInfo.cc
  data/Cluster.cc
)

set(ANALYSIS_INPUT
  input/goat/detail/InputModule.cc
  input/goat/detail//TreeManager.cc
  input/goat/detail/TrackInput.cc
  input/goat/detail/ParticleInput.cc
  input/goat/detail/TriggerInput.cc
  input/goat/detail/DetectorHitInput.cc
  input/goat/detail/TaggerInput.cc
  input/goat/GoatReader.cc
  input/ant/detail/Convert.cc
  input/ant/AntReader.cc
  input/pluto/PlutoReader.cc
  input/pluto/detail/PlutoWrapper.cc
)

set(ANALYSIS_PHYSICS
  physics/test/TestPhysics.h
  physics/test/TestPhysics.cc
  physics/DeltaPlus/DeltaPlusPhysics.h
  physics/DeltaPlus/DeltaPlusPhysics.cc
  physics/omega/omega.h
  physics/omega/omega.cc
  physics/omega/GeoAcceptance.h
  physics/omega/GeoAcceptance.cc
  physics/test/TestAPLCON.h
  physics/test/TestAPLCON.cc
  physics/common/DataOverview.h
  physics/common/DataOverview.cc
  physics/omega/mctrue_acceptance.h
  physics/omega/mctrue_acceptance.cc
  physics/common/CandidatesAnalysis.cc
  )

set(ANALYSIS_ALL
    ${ANALYSIS_BASE}
    ${ANALYSIS_PLOT}
    ${ANALYSIS_DATASTRUCTURES}
    ${ANALYSIS_INPUT}
    ${ANALYSIS_PHYSICS}
)

add_library(analysis ${ANALYSIS_ALL})
target_link_libraries(analysis
  ${ROOT_LIBRARIES}
  ${PLUTO_LIBRARY}
  ${APLCON_LIBRARY}
  base
  tree
  )
