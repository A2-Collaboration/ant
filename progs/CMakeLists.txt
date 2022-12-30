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
add_ant_executable(Ant-completion)

if(AntProgs_MCTools)
    add_ant_executable(Ant-pluto detail/McAction.h)
    add_ant_executable(Ant-mcgun detail/McAction.h)
    add_ant_executable(Ant-cocktail)
    add_ant_executable(Ant-mcdatabase-viewer)
    add_ant_executable(Ant-addTID)
    add_ant_executable(Ant-mc-pi0gun)
    add_ant_executable(Ant-mc-prefilter)
    add_ant_executable(Ant-mc-corrections detail/Corrections.cc)
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
    add_ant_executable(Ant-smoothTaggEff)
    add_ant_executable(Ant-makeLinPol)
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
    add_ant_executable(taggerEnergies)
    add_ant_executable(Ant-repair-keys)
endif()

option(AntProgs_XmasCB "Crystal Ball Christmas ornaments" OFF)

if(AntProgs_XmasCB)
    add_ant_executable(XmasCB)
endif()

option(AntProgs_EtapOmegaG "Build EtaPrime Omega Gamma tools" OFF)

if(AntProgs_EtapOmegaG)
    add_ant_executable(EtapOmegaG_fit)
    target_link_libraries(EtapOmegaG_fit RooFit RooFitCore)
endif()

option(AntProgs_OmegaEtaG "Build Omega Eta Gamma tools" OFF)

if(AntProgs_OmegaEtaG)
    add_ant_executable(OmegaEtaG_fit)
    target_link_libraries(OmegaEtaG_fit RooFit RooFitCore)
endif()

option(AntProgs_Production "Build tools for Production Crosssections" OFF)
if(AntProgs_Production)
    add_ant_executable(SinglePi0_fit)
    add_ant_executable(SigmaPlus_fit detail/SigmaPlus_tools.cc detail/SigmaPlus_tools.h)
    target_link_libraries(SigmaPlus_fit RooFit RooFitCore)
endif()

option(AntProgs_TPCSim "TPC Resolution studies" OFF)
if(AntProgs_TPCSim)
        add_ant_executable(TPCSim detail/TPCSim_tools.cc)
endif()

option(AntProgs_EtapDalitz "Build EtaPrime Dalitz Tools" OFF)

if(AntProgs_EtapDalitz)
    add_ant_executable(EtapDalitz_fit)
    # check if C++17 standard is supported by the compiler
    CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
    if(COMPILER_SUPPORTS_CXX17)
      set_property(TARGET EtapDalitz_fit PROPERTY CXX_STANDARD 17)
      message(STATUS "C++17 will be used for EtapDalitz_fit")
    else()
      message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
    endif()
    # GCC versions prior to the 9.0 release with C++17 support need to be linked to stdc++fs for filesystem support
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
      target_link_libraries(EtapDalitz_fit RooFit RooFitCore roofit_extensions stdc++fs)
    else()
      target_link_libraries(EtapDalitz_fit RooFit RooFitCore roofit_extensions)
    endif()
endif()

option(AntProgs_PCA "Build Principal Component Analysis Tools" OFF)

if(AntProgs_PCA)
    add_ant_executable(PCA)
endif()

option(AntProgs_Pi0TFF "Build Pi0 TFF program" OFF)

if(AntProgs_Pi0TFF)
    add_ant_executable(Pi0TFF)
endif()
