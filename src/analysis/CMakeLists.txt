include_directories(
  .
)


set(ANALYSIS_BASE
  physics/Physics.cc
  utils/combinatorics.h
  utils/matcher.h
  A2GeoAcceptance.cc
  input/DataReader.cc
  utils/ParticleID.cc
  utils/root-addons.cc
  )

set(ANALYSIS_PLOT
  plot/Histogram.cc
  plot/root_draw.cc
  plot/SmartHist.cc
  plot/HistogramFactories.cc
)

set(ANALYSIS_DATASTRUCTURES
  data/Candidate.cc
  data/Particle.cc
  data/Event.cc
  data/TaggerHit.cc
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
  physics/test/TestPhysics.cc
  physics/DeltaPlus/DeltaPlusPhysics.cc
  physics/DeltaTreeGenerator/DeltaTreeGenerator.cc
  physics/omega/omega.cc
  physics/omega/GeoAcceptance.cc
  physics/test/TestAPLCON.cc
  physics/common/DataOverview.cc
  physics/omega/mctrue_acceptance.cc
  physics/common/CandidatesAnalysis.cc
  physics/common/ParticleIDCheck.cc
  physics/common/ProtonCheck.cc
  physics/common/DebugPhysics.cc
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
