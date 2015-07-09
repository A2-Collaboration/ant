# analysis needs some more libraries
find_package(Pluto REQUIRED)
find_package(APLCON REQUIRED)

include_directories(
  .
  ${PLUTO_INCLUDE_DIR}
  ${APLCON_INCLUDE_DIR}
  ../src
  ../src/analysis
  ../src/base
  )

aux_source_directory("." PROGS)
foreach(PROG_SRC ${PROGS})

  get_filename_component(name ${PROG_SRC} NAME_WE)
  message(STATUS "Adding ant executable: ${name}")
  add_executable(${name}
    ${PROG_SRC}
    )

  target_link_libraries(${name} analysis base calibration unpacker)
endforeach()
