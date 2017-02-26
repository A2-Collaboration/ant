include_directories(.)

set(MC_DATANASE
    ProductionDataBase.cc
    Query.cc
)


set(MC_DETAIL
  detail/gp_pPi0.h
)

add_library( mcdatabase
    ${MC_DATANASE}
    ${MC_DETAIL}
)

target_link_libraries(mcdatabase
  MathMore
  base
  tree
)
