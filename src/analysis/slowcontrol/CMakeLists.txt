include_directories(.)

set(SLOWCONTROL
  event_t.h
  SlowControlManager.cc
  SlowControlManager.h
)

set(SLOWCONTROL_PROCESSORS
  SlowControlProcessors.cc
  processors/Processor.h
  processors/AcquScalerProcessor.cc
  processors/ExpTrigger.cc
  processors/Beampolmon.cc
  processors/EPT_Scalers.h
  processors/Tagger_Scalers.h
  processors/PairSpec.cc
)

set(SLOWCONTROL_VARIABLES
  SlowControlVariables.cc
  variables/Variable.h
  variables/TaggerScalers.cc
  variables/Clocks.cc
  variables/Beam.cc
  variables/Pairspec.cc
  variables/Trigger.cc
  variables/TaggEff.cc
)

# since the processors and the variables are
# managed in static registries, they are compiled
# as separate libraries

add_library(slowcontrol_processors  ${SLOWCONTROL_PROCESSORS})
target_link_libraries(slowcontrol_processors expconfig)

add_library(slowcontrol_variables ${SLOWCONTROL_VARIABLES})
target_link_libraries(slowcontrol_variables slowcontrol_processors)

add_library(slowcontrol ${SLOWCONTROL})
target_link_libraries(slowcontrol
    slowcontrol_processors
    slowcontrol_variables
)
