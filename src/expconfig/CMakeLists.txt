set(SRCS
  ExpConfig.cc
  ExpConfig.h
  Detector_t.h
  setups/Setup_2014_EtaPrime.h
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
  detectors/detail/CB_elements.h
  detectors/detail/TAPS_2013_BaF2_elements.h
  detectors/detail/TAPS_2013_PbWO4_elements.h
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

add_library(expconfig ${SRCS})
target_link_libraries(expconfig expconfig-detectors calibration)
