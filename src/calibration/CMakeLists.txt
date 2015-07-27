# analysis needs some more libraries

include_directories(
  .
  ../analysis
  )


set(SRCS
    CalibrationManager.h
    Calibration.cc
    modules/Timing.cc
    modules/Integral.cc
    modules/Scaler.cc
    modules/EnergyInvariantMass.cc
    modules/TAPS_ShowerCorrection.h
    modules/CB_Energy.cc
    modules/TAPS_Energy.cc
    converters/GeSiCa_SADC.h
    converters/CATCH_TDC.cc
    converters/MultiHit16bit.h
    converters/ScalerFrequency.cc
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration base analysis)
