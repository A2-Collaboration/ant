# analysis needs some more libraries

include_directories(
  .
  ../analysis
  )


set(SRCS
    Calibration.h
    Editor.cc
    DataBase.cc
    DataManager.cc
    detail/Helpers.h
    modules/Time.cc
    modules/Energy.cc
    modules/Scaler.cc
    modules/TAPS_ShowerCorrection.h
    modules/CB_Energy.cc
    modules/TAPS_Energy.cc
    modules/PID_Energy.cc
    modules/PID_PhiAngle.cc
    modules/TAPSVeto_Energy.cc
    converters/GeSiCa_SADC.h
    converters/CATCH_TDC.cc
    converters/MultiHit16bit.h
    converters/ScalerFrequency.cc
    gui/Manager.cc
    gui/Manager_traits.h
    gui/BaseFunctions.cc
    gui/GUIbase.h
    gui/FitFunction.cc
    gui/CalCanvas.cc
    gui/GUIElements.cc
    gui/TF1Knobs.cc
    gui/FitGausPol3.cc
    gui/AvgBuffer.h
    gui/ACECanvas.cc
    gui/Dialogs.cc
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration base analysis expconfig-detectors)
