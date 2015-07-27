# analysis needs some more libraries

include_directories(
  .
  ../analysis
  )


set(SRCS
    CalibrationManager.h
    Calibration.h
    Calibration.cc
    modules/Timing.h
    modules/Timing.cc
    modules/Integral.h
    modules/Integral.cc
    modules/Scaler.h
    modules/Scaler.cc
    modules/EnergyInvariantMass.h
    modules/EnergyInvariantMass.cc
    modules/TAPS_ShowerCorrection.h
    converters/GeSiCa_SADC.h
    converters/CATCH_TDC.h
    converters/CATCH_TDC.cc
    converters/MultiHit16bit.h
    converters/ScalerFrequency.h
    converters/ScalerFrequency.cc    
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration base analysis)
