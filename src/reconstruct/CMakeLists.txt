set(SRCS
  Reconstruct.h
  Reconstruct.cc
  Reconstruct_traits.h
  HitMatching.h
  ApplyCalibrations.h
  Clustering.h
  )

include_directories(.)
add_library(reconstruct ${SRCS})
target_link_libraries(reconstruct base unpacker calibration)
