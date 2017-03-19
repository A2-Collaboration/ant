include_directories(.)
include_directories(../analysis)

set(SRCS
    Calibration.h
    DataBase.cc
    DataManager.cc
    Editor.cc
    modules/Time.cc
    modules/Energy.cc
    modules/CalibType.cc
    modules/Energy_GUI.cc
    modules/TAPS_ShowerCorrection.h
    modules/CB_Energy.cc
    modules/CB_TimeWalk.cc
    modules/TAPS_Energy.cc
    modules/TAPS_ShortEnergy.cc
    modules/PID_Energy.cc
    modules/PID_PhiAngle.cc
    modules/TAPSVeto_Energy.cc
    modules/TAPS_Time.cc
    modules/TAPS_ToF.cc
    modules/TAPSVeto_Time.cc
    modules/TaggEff.cc
    modules/CB_SourceCalib.cc
    modules/Tagger_QDC.cc
    modules/ClusterCorrection.cc
    converters/GeSiCa_SADC.h
    converters/MultiHit.h
    converters/MultiHitReference.cc
    converters/CATCH_TDC.h
    gui/Manager.cc
    gui/Manager_traits.h
    gui/ManagerWindow.cc
    gui/ManagerWindow_traits.h
    gui/EditorWindow.cc
    gui/EditorCanvas.cc
    gui/Indicator_traits.h
    gui/Indicator.cc
    gui/CalCanvas.cc
    gui/AvgBuffer.h
    gui/AvgBuffer_traits.h
    gui/Dialogs.cc
    fitfunctions/BaseFunctions.cc
    fitfunctions/KnobsTF1.cc
    fitfunctions/FitFunction.cc
    fitfunctions/FitGaus.cc
    fitfunctions/FitGausPol0.cc
    fitfunctions/FitGausPol1.cc
    fitfunctions/FitGausPol3.cc
    fitfunctions/FitTimewalk.cc
    fitfunctions/FitLandau.cc
    fitfunctions/FitLandauPol0.cc
    fitfunctions/FitLandauExpo.cc
    fitfunctions/FitGausexpo.cc
    fitfunctions/FitWeibullLandauPol1.cc
    modules/detail/TH2Storage.cc
  )


add_library(calibration ${SRCS})
target_link_libraries(calibration
    base tree expconfig_detectors cbtaps_display
    Spectrum MathMore
    )
