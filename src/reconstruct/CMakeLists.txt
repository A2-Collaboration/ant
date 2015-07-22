set(SRCS
  Reconstruct.h
  Reconstruct.cc
  Reconstruct_traits.h
  Clustering.h
  TrackBuilder.h
  TrackBuilder.cc
  )

include_directories(.)
add_library(reconstruct ${SRCS})
target_link_libraries(reconstruct base unpacker)
