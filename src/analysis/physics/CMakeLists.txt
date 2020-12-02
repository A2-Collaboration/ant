set(LIBS "")

macro(add_ant_physics SRCFILE)
  string(REGEX REPLACE "\\.cc$" "" BASENAME ${SRCFILE})
  string(REGEX REPLACE "[^A-Za-z0-9_]" "_" BASENAME ${BASENAME})
  set(LIBNAME "physicsclass_${BASENAME}")
  add_library(${LIBNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${SRCFILE})
  target_link_libraries(${LIBNAME} third_party_interface)
  # the additional linker flag is very important if
  # physics classes define symbols (such as plotter classes)
  # with the same name (see for example tunings/TriggerSimulation and check/PIDEfficiencyCheck)
  # the linker then picks just one, and this leads to weird segfaults in Ant-plot
  set_target_properties(${LIBNAME} PROPERTIES
    LINK_FLAGS "-Wl,-Bsymbolic"
    # tweak also the output into separare lib dir
    OUTPUT_NAME ${BASENAME}
    PREFIX "" # no lib prefix
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/PhysicsClasses"
  )
  list(APPEND LIBS ${LIBNAME})
  set(LIBS ${LIBS} PARENT_SCOPE)
endmacro()

add_subdirectory(calibration)
add_subdirectory(check)
add_subdirectory(common)
add_subdirectory(etaprime)
add_subdirectory(misc)
add_subdirectory(omega)
add_subdirectory(pi0)
add_subdirectory(production)
add_subdirectory(scratch)
add_subdirectory(test)
add_subdirectory(tpc)
add_subdirectory(tunings)
add_subdirectory(Compton)

add_library(physics_classes INTERFACE)
target_link_libraries(physics_classes INTERFACE ${LIBS})
