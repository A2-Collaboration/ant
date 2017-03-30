include_directories(
  ${CMAKE_SOURCE_DIR}/src
  )

macro(add_ant_executable name)
  set(EXEFILE "${name}.cc")
  add_executable(${name} ${EXEFILE} ${ARGN})
  target_link_libraries(${name} analysis base calibration unpacker reconstruct pluto)
endmacro()

option(AntProgs_CalibTools  "Calibration Tools" ON)
option(AntProgs_MCTools     "Monte Carlo Tools" ON)
option(AntProgs_TuningTools "Tuning Tools"      ON)
option(AntProgs_DebugTools  "Debug Tools"       ON)

add_ant_executable(Ant)
add_ant_executable(Ant-plot)

add_ant_executable(Ant-chain)
add_ant_executable(Ant-hadd)
add_ant_executable(Ant-info)

if(AntProgs_MCTools)
    add_ant_executable(Ant-pluto detail/McAction.h)
    add_ant_executable(Ant-mcgun detail/McAction.h)
    add_ant_executable(Ant-cocktail)
    add_ant_executable(Ant-mcdatabase-viewer)
    add_ant_executable(Ant-addTID)
endif()

if(AntProgs_CalibTools)
    add_ant_executable(Ant-calib)
    add_ant_executable(Ant-calib-regedit)
    add_ant_executable(Ant-calib-editor)
    add_ant_executable(Ant-calib-viewer)
    add_ant_executable(Ant-calib-dump)
    add_ant_executable(Ant-calib-readin)
    add_ant_executable(Ant-calib-smooth)
    add_ant_executable(Ant-altVetoCalTool)
    add_ant_executable(Ant-makeTaggEff detail/taggEffClasses.cc)
endif()

if(AntProgs_TuningTools)
    add_ant_executable(Ant-mcsmearing)
    add_ant_executable(Ant-ECorr)
    add_ant_executable(Ant-tuning-viewer)
    add_ant_executable(Ant-makeSigmas)
endif()

if(AntProgs_DebugTools)
    add_ant_executable(Ant-rawdump)
    add_ant_executable(Ant-fakeRaw)
    add_ant_executable(Ant-treeTool)
    add_ant_executable(Ant-copyTree)
    add_ant_executable(Ant-compareHists)
endif()

option(AntProgs_SimpleTools "Simple Tools" OFF)

if(AntProgs_SimpleTools)
    add_ant_executable(processEtapTree)
    add_ant_executable(compare_ant_goat)
    add_ant_executable(compare_tree_events)
    add_ant_executable(simulate_livetime)
    add_ant_executable(taggEffPlot detail/taggEffClasses.cc detail/tools.h)
endif()

option(AntProgs_XmasCB "Crystal Ball Christmas ornaments" OFF)

if(AntProgs_XmasCB)
    add_ant_executable(XmasCB)
endif()

option(AntProgs_EtapOmegaG "Build EtaPrime Omega Gamma tools" OFF)

if(AntProgs_EtapOmegaG)
    add_ant_executable(EtapOmegaG_fit)
endif()
