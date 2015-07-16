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
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration
  base
  analysis
  )
