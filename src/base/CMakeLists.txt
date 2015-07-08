include_directories(.)

set(SRCS
  types.cc
  types.h
  interval.h
  printable.h
  printable.cc
  Logger.cc
  Logger.h
  detail/easylogging++.h
  tmpfile_t.h
  tmpfile_t.cc
  stl_helpers.h
)

add_library(base ${SRCS})
