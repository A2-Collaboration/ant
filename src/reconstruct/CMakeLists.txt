set(SRCS
  Reconstruct.h
  Reconstruct.cc
  Reconstruct_traits.h
  Clustering.h
  Clustering.cc
  CandidateBuilder.h
  CandidateBuilder.cc
  UpdateableManager.h
  UpdateableManager.cc
  detail/Clustering_NextGen.h
  )

include_directories(.)
add_library(reconstruct ${SRCS})
target_link_libraries(reconstruct base unpacker)
