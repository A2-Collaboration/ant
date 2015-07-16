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
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration
  base
  analysis
  )
