set(SRCS
  ExpConfig.cc
  ExpConfig.h
  detectors/TAPS.h
  detectors/CB.h
  detectors/CB.cc
  setups/Setup_2014_EtaPrime.h
  )

include_directories(.)
add_library(expconfig ${SRCS})
target_link_libraries(expconfig base)
