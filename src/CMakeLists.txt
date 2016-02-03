# this makes cross-includes easy,
# for example for like
# >>  #include "base/Logger.h"
# but don't forget to link against
# the corresponding library then
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# for generated files by configure_file()
include_directories(${CMAKE_CURRENT_BINARY_DIR})


macro(add_ant_subdirectory name)
  message(STATUS "Configuring ${name}...")
  add_subdirectory(${name})
endmacro()

macro(append_header_files varname)
  set(more_headers "")
  foreach(f ${${varname}})
    get_filename_component(ext ${f} EXT)
    if(${ext} STREQUAL ".cc")
      get_filename_component(dir ${f} PATH)
      get_filename_component(name_we ${f} NAME_WE)
      set(f_header "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name_we}.h")
      if(EXISTS ${f_header})
        list(APPEND more_headers ${f_header})
      endif()
    endif()
  endforeach()
  list(LENGTH more_headers n_headers)
  if(${n_headers} GREATER 0)
    message(STATUS "Added ${n_headers} extra header files")
  endif()  
  list(APPEND ${varname} ${more_headers})
endmacro()

# the order is important if libraries from
# other subdirectories are needed
add_ant_subdirectory(base)
add_ant_subdirectory(expconfig)
add_ant_subdirectory(tree)
add_ant_subdirectory(unpacker)
add_ant_subdirectory(reconstruct)
add_ant_subdirectory(analysis)
add_ant_subdirectory(calibration)
add_ant_subdirectory(simulation/mc)
add_ant_subdirectory(root-addons)
