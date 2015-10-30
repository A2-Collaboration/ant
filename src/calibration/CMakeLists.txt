include_directories(.)
include_directories(../analysis)

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
    modules/TAPS_ShortEnergy.cc
    modules/PID_Energy.cc
    modules/PID_PhiAngle.cc
    modules/TAPSVeto_Energy.cc
    modules/TAPS_Time.cc
    converters/GeSiCa_SADC.h
    converters/CATCH_TDC.cc
    converters/MultiHit16bit.h
    converters/ScalerFrequency.cc
    gui/Manager.cc
    gui/Manager_traits.h
    gui/ManagerWindow.cc
    gui/EditorCanvas.cc
    gui/EditorWindow.cc
    gui/Indicator_traits.h
    gui/Indicator.cc
    gui/CalCanvas.cc
    gui/AvgBuffer.h
    gui/Dialogs.cc
    fitfunctions/BaseFunctions.cc
    fitfunctions/KnobsTF1.cc
    fitfunctions/FitFunction.cc
    fitfunctions/FitGaus.cc
    fitfunctions/FitGausPol0.cc
    fitfunctions/FitGausPol3.cc
    fitfunctions/FitTimewalk.cc
  )

add_library(calibration ${SRCS})
target_link_libraries(calibration base analysis expconfig-detectors)
