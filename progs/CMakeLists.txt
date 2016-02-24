include_directories(
  ${CMAKE_SOURCE_DIR}/src
  )

macro(add_ant_executable name)
  set(EXEFILE "${name}.cc")
  add_executable(${name} ${EXEFILE})
  target_link_libraries(${name} analysis base calibration unpacker reconstruct mc)
endmacro()

add_ant_executable(Ant)
add_ant_executable(Ant-pluto)
add_ant_executable(Ant-cocktail)

add_ant_executable(Ant-calib)
add_ant_executable(Ant-calib-regedit)
add_ant_executable(Ant-calib-editor)
add_ant_executable(Ant-calib-viewer)

add_ant_executable(Ant-chain)
add_ant_executable(Ant-addTID)
add_ant_executable(Ant-rawdump)
add_ant_executable(Ant-info)
add_ant_executable(Ant-treeTool)

option(OmegaTools "Build Omega Tools" OFF)

if(OmegaTools)
endif()

option(SimpleTools "Simple Tools" OFF)

if(SimpleTools)
    add_ant_executable(show_cocktail_database)
    add_ant_executable(compare_ant_goat)
endif()

option(XmasCB "Crystal Ball Christmas ornaments" OFF)

if(XmasCB)
    add_ant_executable(XmasCB)
endif()

option(EtapOmegaG "Build EtaPrime Omega Gamma tools" OFF)

if(EtapOmegaG)
    add_ant_executable(EtapOmegaG_plot)
endif()
