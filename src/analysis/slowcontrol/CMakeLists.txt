include_directories(.)

set(SLOWCONTROL
  SlowControl.cc
  SlowControlCreator.cc
  SlowControlManager.cc
)

set(SLOWCONTROL_PROCESSORS
  SlowControlProcessors.h
  SlowControlProcessors.cc
  processors/Processor.h
  processors/EPT_Scalers.h
  processors/TotalLivetime.h
  processors/Tagger_Scalers.h
  processors/Tagger_Livetime.h
  processors/PairSpec.h
)

set(SLOWCONTROL_VARIABLES
  SlowControlVariables.cc
  variables/Variable.h
  variables/TaggerScalers.cc
  variables/TotalLivetime.h
  variables/BeamCurrent.h
  variables/PhotonFlux.h
)

# since the processors and the variables are
# managed in static registries, they are compiled
# as separate libraries

add_library(slowcontrol_processors  ${SLOWCONTROL_PROCESSORS})

append_header_files(SLOWCONTROL_VARIABLES)
add_library(slowcontrol_variables ${SLOWCONTROL_VARIABLES})


add_library(slowcontrol ${SLOWCONTROL})
target_link_libraries(slowcontrol
    slowcontrol_processors
    slowcontrol_variables
)