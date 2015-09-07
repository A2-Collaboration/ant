include_directories(
	.
	${CMAKE_CURRENT_BINARY_DIR}
)

#set data file, parser and output for precompiled cocktail-database
set(XSECTION_DATA      ${CMAKE_CURRENT_SOURCE_DIR}/xsections.dat)
set(XSECTION_PARSER    ${CMAKE_CURRENT_SOURCE_DIR}/DecayParser.py)
set(XSECTION_DATA_SRC  ${CMAKE_CURRENT_BINARY_DIR}/chgen.h)

set(SIMULATION_MC
  A2Channels.cc
  A2Cocktail.cc
)


# use Python-Script for precompiled database
ADD_CUSTOM_COMMAND(
	OUTPUT ${XSECTION_DATA_SRC}
	COMMAND ${XSECTION_PARSER} ${XSECTION_DATA} ${XSECTION_DATA_SRC}
	DEPENDS ${XSECTION_DATA} ${XSECTION_PARSER}
)

add_library( mc 
  ${SIMULATION_MC}
  ${XSECTION_DATA_SRC}
)

target_link_libraries( mc
  ${ROOT_LIBRARIES}
  ${PLUTO_LIBRARY}
  MathMore
)
