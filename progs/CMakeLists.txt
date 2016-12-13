include_directories(
  ${CMAKE_SOURCE_DIR}/src
  )

macro(add_ant_executable name)
  set(EXEFILE "${name}.cc")
  add_executable(${name} ${EXEFILE} ${ARGN})
  target_link_libraries(${name} analysis base calibration unpacker reconstruct pluto)
endmacro()

add_ant_executable(Ant)
add_ant_executable(Ant-pluto)
add_ant_executable(Ant-cocktail)

add_ant_executable(Ant-calib)
add_ant_executable(Ant-calib-regedit)
add_ant_executable(Ant-calib-editor)
add_ant_executable(Ant-calib-viewer)

add_ant_executable(Ant-mcsmearing)

add_ant_executable(Ant-makeTaggEff detail/taggEffClasses.cc)

add_ant_executable(Ant-chain)
add_ant_executable(Ant-hadd)
add_ant_executable(Ant-addTID)
add_ant_executable(Ant-rawdump)
add_ant_executable(Ant-fakeRaw)
add_ant_executable(Ant-info)
add_ant_executable(Ant-treeTool)
add_ant_executable(Ant-copyTree)
add_ant_executable(Ant-makeSigmas)
add_ant_executable(Ant-compareHists)
add_ant_executable(Ant-altVetoCalTool)

option(AntProgs_Omega "Build Omega Analysis Tools" OFF)

if(AntProgs_Omega)
    add_ant_executable(OmegaEtaG2_plot)
endif()

option(AntProgs_Prod "Build Production Analysis Tools" OFF)
if (AntProgs_Prod)
    add_ant_executable(TriplePi0_plot)
    add_ant_executable(TriplePi0_cuts)
endif()

option(AntProgs_SimpleTools "Simple Tools" OFF)

if(AntProgs_SimpleTools)
    add_ant_executable(processEtapTree)
    add_ant_executable(show_cocktail_database)
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
    add_ant_executable(EtapOmegaG_plot)
    add_ant_executable(EtapOmegaG_fit)
endif()

option(AntProgs_Dalitz "Build eta/eta' Dalitz analysis tools" OFF)

if(AntProgs_Dalitz)
    add_ant_executable(EtaDalitz_plot)
    add_ant_executable(EtapDalitz_plot)
endif()
