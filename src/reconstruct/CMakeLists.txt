include_directories(.)

set(SRCS
  Reconstruct.cc
  Reconstruct_traits.h
  Clustering.cc
  CandidateBuilder.cc
  UpdateableManager.cc
  detail/Clustering_NextGen.h
  )

append_header_files(SRCS)

add_library(reconstruct ${SRCS})
target_link_libraries(reconstruct base expconfig)
