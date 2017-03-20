set(SRCS
  event_t.cc
  DataReader.h
  goat/GoatReader.cc
  ant/AntReader.cc
  pluto/PlutoReader.cc
  pluto/detail/PlutoWrapper.cc
)

add_library(analysis_input ${SRCS})
target_link_libraries(analysis_input third_party_interface)
