# analysis needs some more libraries
include_directories(
  ../src
  )

aux_source_directory("." PROGS)
foreach(PROG_SRC ${PROGS})

  get_filename_component(name ${PROG_SRC} NAME_WE)
  message(STATUS "Adding ant executable: ${name}")
  add_executable(${name}
    ${PROG_SRC}
    )

  target_link_libraries(${name} analysis base calibration unpacker reconstruct)
endforeach()
