set(SRCS
  ExpConfig.cc
  ExpConfig.h
  Detector_t.h
  detectors/TAPS.h
  detectors/TAPS.cc
  detectors/CB.h
  detectors/CB.cc
  detectors/detail/CB_elements.h
  setups/Setup_2014_EtaPrime.h
  )

include_directories(. ../analysis)
add_library(expconfig ${SRCS})
target_link_libraries(expconfig base tree)
