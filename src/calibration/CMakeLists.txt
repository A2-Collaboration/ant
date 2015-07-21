# analysis needs some more libraries

include_directories(
  .
  ../analysis
  )


set(SRCS
    Calibration.h
    Calibration.cc
    modules/EnergyInvariantMass.h
    modules/EnergyInvariantMass.cc
    modules/TimingCATCH.cc
    modules/TimingCATCH.h
    modules/IntegralSADC.cc
    modules/IntegralSADC.h
    modules/IntegralCAEN.cc
    modules/IntegralCAEN.h
    modules/TimingTAPS.cc
    modules/TimingTAPS.h
    modules/IntegralTAPS.cc
    modules/IntegralTAPS.h
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration
  base
  analysis
  )
