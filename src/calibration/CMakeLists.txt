# analysis needs some more libraries

include_directories(
  .
  ../analysis
  )


set(CALIBRATION_BASE
    CalibrationApply.h
    BaseCalModule.h
    BaseCalModule.cc
    TestCalCB.h
    TestCalCB.cc
  )

add_library(calibration ${CALIBRATION_BASE})
target_link_libraries(calibration
  ${ROOT_LIBRARIES}
  base
  analysis
  )

