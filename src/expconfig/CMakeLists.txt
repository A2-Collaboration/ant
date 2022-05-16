include_directories(.)
include_directories(../analysis ../calibration)

set(SRCS
  ExpConfig.cc
  )

set(SRCS_SETUPS
  setups/Setup.cc
  setups/SetupRegistry.cc
  setups/Setup_2012_12_Compton.cc
  setups/Setup_2014_EPT.cc
  setups/Setup_2014_07_EPT_Prod.cc
  setups/Setup_2014_10_EPT_Prod.cc
  setups/Setup_2014_12_EPT_Prod.cc
  setups/Setup_2015_01_Pion.cc
  setups/Setup_CBSourceCalib.cc
  setups/Setup_2007_Base.cc
  setups/Setup_2007_06.cc
  setups/Setup_2007_07.cc
  setups/Setup_2017_03.cc
  setups/Setup_2017_05.cc
  setups/Setup_2017Plus_NewTagger_Base.cc
  setups/Setup_2017_12.cc
  setups/Setup_2018_03.cc
  setups/Setup_2018_05.cc
  setups/Setup_2018_09.cc
  setups/Setup_2018_11.cc
  setups/Setup_2019_01.cc
  setups/Setup_2019_06.cc
  setups/Setup_2019_07.cc
  setups/Setup_2019_09.cc
  setups/Setup_2021_10.cc
  setups/Setup_2022_05.cc
  setups/Setup_2010_03_Base.cc
  setups/Setup_2010_09_Compton.cc
  )

set(SRCS_DETECTORS
  detectors/Trigger.cc
  detectors/EPT.cc
  detectors/Tagger.cc
  detectors/TAPS.cc
  detectors/CB.cc
  detectors/PID.cc
  detectors/TAPSVeto.cc
  detectors/Cherenkov.h
  detectors/detail/EPT_2014_elements.h
  detectors/detail/Tagger_2007_1508_elements.h
  detectors/detail/Tagger_2010_03_450_elements.h
  detectors/detail/Tagger_2015_450_elements.h
  detectors/detail/Tagger_2016_06_1557_elements.h
  detectors/detail/Tagger_2017_12_883_elements.h
  detectors/detail/Tagger_2018_03_883_elements.h
  detectors/detail/Tagger_2018_11_883_elements.h
  detectors/detail/Tagger_2019_01_883_elements.h
  detectors/detail/Tagger_2019_06_450_elements.h
  detectors/detail/Tagger_2019_07_1604_elements.h
  detectors/detail/CB_elements.h
  detectors/detail/TAPS_2007_BaF2_elements.h
  detectors/detail/TAPS_2009_03_BaF2_elements.h
  detectors/detail/TAPS_2009_03_PbWO4_elements.h
  detectors/detail/TAPS_2013_11_BaF2_elements.h
  detectors/detail/TAPS_2013_11_PbWO4_elements.h
  detectors/detail/TAPS_2021_10_BaF2_elements.h
  detectors/detail/TAPS_2021_10_PbWO4_elements.h
  detectors/detail/TAPSVeto_2007_BaF2_elements.h
  detectors/detail/TAPSVeto_2009_03_BaF2_elements.h
  detectors/detail/TAPSVeto_2009_03_PbWO4_elements.h
  detectors/detail/TAPSVeto_2013_11_BaF2_elements.h
  detectors/detail/TAPSVeto_2013_11_PbWO4_elements.h
  detectors/detail/TAPSVeto_2014_BaF2_elements.h
  detectors/detail/TAPSVeto_2014_PbWO4_elements.h
  detectors/detail/PID_2004_elements.h
  detectors/detail/PID_2007_elements.h
  detectors/detail/PID_2009_05_elements.h
  detectors/detail/PID_2009_06_elements.h
  detectors/detail/PID_2009_07_elements.h
  detectors/detail/PID_2014_elements.h
)


# make detector stuff own library,
# very important to get init of static members at correct runtime
add_library(expconfig_detectors ${SRCS_DETECTORS})
target_link_libraries(expconfig_detectors base ${ROOT_LIBRARIES})

add_library(expconfig_setups ${SRCS_SETUPS})
target_link_libraries(expconfig_setups expconfig_detectors calibration)

add_library(expconfig ${SRCS})
target_link_libraries(expconfig expconfig_setups)
