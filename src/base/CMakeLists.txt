set(ANT_PATH_DATABASE "${CMAKE_SOURCE_DIR}/database")
set(ANT_PATH_GITREPO  "${CMAKE_SOURCE_DIR}")

configure_file(Paths.h.in Paths.h)

set(SRCS_BASE
  Paths.h
  types.h
  Detector_t.cc
  ParticleType.cc
  ParticleTypeTree.cc
  interval.h
  piecewise_interval.h
  interval_algo.h
  Logger.cc
  tmpfile_t.cc
  Format.h
  Tree.h
  WrapTFile.cc
  matrixstack.cc
  GitInfo.cc
  OptionsList.cc
  ProgressCounter.cc
  TF1Ext.h
  PlotExt.cc
  WrapTTree.cc
  Interpolator.cc
  Array2D.cc
  TH_ext.cc
  BinSettings.cc
  ClippedInterpolatorWrapper.cc
  FloodFillAverages.h
  bitflag.h
  SavitzkyGolay.cc
  PhysicsMath.h
  )

set(SRCS_VEC
  vec/vec3.cc
  vec/vec2.cc
  vec/LorentzVec.cc
)

set(SRCS_MATH_FUNCTIONS
  math_functions/voigtian.h
  math_functions/voigtian.cc
  math_functions/AsymGaus.cc
  math_functions/CrystalBall.cc
  math_functions/Linear.h
)

set(SRCS_STD_EXT
  std_ext/math.cc
  std_ext/memory.h
  std_ext/misc.h
  std_ext/string.h
  std_ext/time.h
  std_ext/container.h
  std_ext/system.cc
  std_ext/convert.h
  std_ext/iterators.h
  std_ext/mapped_vectors.h
  std_ext/shared_ptr_container.h
  std_ext/printable.h
  std_ext/variadic.h
  std_ext/vector.h
)

set(SRCS
  ${SRCS_BASE}
  ${SRCS_MATH_FUNCTIONS}
  ${SRCS_STD_EXT}
  ${SRCS_VEC}
)

add_library(base ${SRCS})
target_link_libraries(base third_party ${ROOT_LIBRARIES} ${GSL_LIBRARIES} ${PLUTO_LIBRARY})
