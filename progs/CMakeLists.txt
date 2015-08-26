include_directories(
  ../src
  )

macro(add_ant_executable name)
  set(EXEFILE "${name}.cc")
  add_executable(${name} ${EXEFILE})
  target_link_libraries(${name} analysis base calibration unpacker reconstruct)
endmacro()

add_ant_executable(Ant)
add_ant_executable(Ant-calib)
add_ant_executable(Ant-calib-editor)
