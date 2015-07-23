set(SRCS
  ExpConfig.cc
  ExpConfig.h
  Detector_t.h
  )

set(SRCS_SETUPS
  Setup.h
  Setup.cc
  setups/Setup_2014_EtaPrime.cc
  setups/Setup_2015_dummy.cc
  )

set(SRCS_DETECTORS
  detectors/Trigger.h
  detectors/EPT.h
  detectors/EPT.cc
  detectors/TAPS.h
  detectors/TAPS.cc
  detectors/CB.h
  detectors/CB.cc
  detectors/PID.h
  detectors/PID.cc
  detectors/TAPSVeto.h
  detectors/TAPSVeto.cc
  detectors/detail/EPT_2014_elements.h
  detectors/detail/CB_elements.h
  detectors/detail/TAPS_2013_BaF2_elements.h
  detectors/detail/TAPS_2013_PbWO4_elements.h
  detectors/detail/TAPSVeto_2013_PbWO4_elements.h
  detectors/detail/TAPSVeto_2013_BaF2_elements.h
  detectors/detail/TAPSVeto_2014_BaF2_elements.h
  detectors/detail/PID_2004_elements.h
  detectors/detail/PID_2009_05_elements.h
  detectors/detail/PID_2009_06_elements.h
  detectors/detail/PID_2009_07_elements.h
  detectors/detail/PID_2014_elements.h
)

include_directories(. ../analysis ../calibration)

# make detector stuff own library,
# very important to get init of static members at correct runtime
add_library(expconfig-detectors ${SRCS_DETECTORS})
target_link_libraries(expconfig-detectors ${ROOT_LIBRARIES})

add_library(expconfig-setups ${SRCS_SETUPS})
target_link_libraries(expconfig-setups expconfig-detectors calibration)

add_library(expconfig ${SRCS})
target_link_libraries(expconfig expconfig-setups)
