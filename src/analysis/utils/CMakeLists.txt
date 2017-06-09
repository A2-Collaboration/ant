set(SRCS
  Combinatorics.h
  Matcher.h
  A2GeoAcceptance.cc
  ParticleID.cc
  RootAddons.cc
  ParticleTools.cc
  TimeSmearingHack.cc
  fitter/Fitter.cc
  fitter/KinFitter.cc
  fitter/TreeFitter.cc
  Uncertainties.cc
  ProtonPermutation.cc
  ClusterTools.cc
  MCSmear.cc
  MCFakeReconstructed.cc
  PullsWriter.h
  uncertainties/Constant.cc
  uncertainties/MCExtracted.cc
  uncertainties/Optimized.cc
  uncertainties/Interpolated.cc
  uncertainties/FitterSergey.cc
  uncertainties/MeasuredProton.cc
  uncertainties/MCSmearingAdlarson.cc
  MCWeighting.cc
  TriggerSimulation.cc
  ProtonPhotonCombs.cc
  )

add_library(analysis_utils ${SRCS})
target_link_libraries(analysis_utils ${APLCONpp_LIBRARY} third_party_interface)
