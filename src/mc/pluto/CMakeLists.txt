set(SRCS
  PlutoGenerator.cc
  PlutoExtensions.h
  PlutoFactory.cc
  utils/PlutoTID.cc
)

add_library( pluto  ${SRCS})

target_link_libraries(pluto
  ${ROOT_LIBRARIES}
  ${PLUTO_LIBRARY}
  MathMore
  base
  tree
  mcdatabase
)
