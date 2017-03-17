include_directories(.)

set(MC_DATABASE
    ProductionDataBase.cc
    Query.cc
)


set(MC_DETAIL
  detail/gp_pPi0.h
  detail/gp_pg.h
  detail/gp_nPiP.h
  detail/gp_pEta.h
  detail/gp_pEtaPi0.h
  detail/gp_pEtaPrime.h
  detail/gp_pOmega.h
  detail/gp_pPi0Pi0.h
  detail/gp_pRho.h
  detail/gp_SigmaPlusK0S.h
)

add_library( mcdatabase
    ${MC_DATABASE}
    ${MC_DETAIL}
)

target_link_libraries(mcdatabase
  MathMore
  base
  tree
)
