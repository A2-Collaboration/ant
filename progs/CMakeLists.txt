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

# temporary executables...

add_executable(test_gui test_gui.cc)
target_link_libraries(test_gui ${ROOT_LIBRARIES})

add_executable(test_display test_display.cc)
target_link_libraries(test_display base)