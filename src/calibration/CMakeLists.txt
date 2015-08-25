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
    modules/Time.cc
    modules/Energy.cc
    modules/Scaler.cc
    modules/TAPS_ShowerCorrection.h
    modules/CB_Energy.cc
    modules/CB_TimeWalk.cc
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
    gui/ManagerWindow.cc
    gui/GUIbase.h
    gui/CalCanvas.cc
    gui/GUIElements.cc
    gui/AvgBuffer.h
    gui/ACECanvas.cc
    gui/Dialogs.cc
    fitfunctions/BaseFunctions.cc
    fitfunctions/TF1Knobs.cc
    fitfunctions/FitFunction.cc
    fitfunctions/FitGaus.cc
    fitfunctions/FitGausPol0.cc
    fitfunctions/FitGausPol3.cc
    fitfunctions/FitTimewalk.cc
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration base analysis expconfig-detectors)
