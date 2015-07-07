set(SRCS
  ExpConfig.cc
  ExpConfig.h
  )

add_library(expconfig ${SRCS})
target_link_libraries(expconfig base)
