include_directories(
  ../src
  )

macro(add_ant_executable name)
  set(EXEFILE "${name}.cc")
  add_executable(${name} ${EXEFILE})
  target_link_libraries(${name} analysis base calibration unpacker reconstruct mc)
endmacro()

add_ant_executable(Ant)
add_ant_executable(Ant-calib)
add_ant_executable(Ant-calib-editor)
add_ant_executable(Ant-cocktail)
add_ant_executable(Ant-pluto)
add_ant_executable(Ant-chain)
add_ant_executable(Ant-addTID)
add_ant_executable(Ant-info)

option(OmegaTools "Build Omega Tools" OFF)

if(OmegaTools)
    add_ant_executable(omega_histograms)
    add_ant_executable(omega_tree)
endif()

option(SimpleTools "Simple Tools" OFF)

if(SimpleTools)
    add_ant_executable(show_cocktail_database)
endif()
