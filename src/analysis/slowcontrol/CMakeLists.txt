include_directories(.)

set(SLOWCONTROL
  SlowControlManager.cc
  SlowControlManager.h
)

set(SLOWCONTROL_PROCESSORS
  SlowControlProcessors.cc
  processors/Processor.h
  processors/AcquScalerProcessor.cc
  processors/ExpTrigger.h
  processors/Beampolmon.h
  processors/EPT_Scalers.h
  processors/Tagger_Scalers.h
  processors/PairSpec.h
)

set(SLOWCONTROL_VARIABLES
  SlowControlVariables.cc
  variables/Variable.h
  variables/TaggerScalers.cc
  variables/TotalLivetime.h
  variables/BeamCurrent.h
  variables/PhotonFlux.cc
)

# since the processors and the variables are
# managed in static registries, they are compiled
# as separate libraries

append_header_files(SLOWCONTROL_PROCESSORS)
add_library(slowcontrol_processors  ${SLOWCONTROL_PROCESSORS})
target_link_libraries(slowcontrol_processors expconfig)

append_header_files(SLOWCONTROL_VARIABLES)
add_library(slowcontrol_variables ${SLOWCONTROL_VARIABLES})
target_link_libraries(slowcontrol_variables slowcontrol_processors)

add_library(slowcontrol ${SLOWCONTROL})
target_link_libraries(slowcontrol
    slowcontrol_processors
    slowcontrol_variables
)
