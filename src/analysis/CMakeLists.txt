include_directories(
  .
)

set(SRCS
  physics/Physics.cc
  physics/PhysicsManager.cc
  physics/manager_t.h
  physics/Plotter.cc
  )

add_subdirectory(input)
add_subdirectory(utils)
add_subdirectory(plot)
add_subdirectory(slowcontrol)
add_subdirectory(physics)

add_library(analysis ${SRCS})
target_link_libraries(analysis
  physics_classes
  analysis_utils
  analysis_input
  analysis_plot
  ${ROOT_LIBRARIES}
  HistPainter
  ${PLUTO_LIBRARY}
  ${APLCONpp_LIBRARY}
  base
  tree
  expconfig
  cbtaps_display
  analysis_codes
  slowcontrol
  )
