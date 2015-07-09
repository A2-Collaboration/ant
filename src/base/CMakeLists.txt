include_directories(.)

set(SRCS
  types.cc
  types.h
  interval.h
  piecewise_interval.h
  printable.h
  printable.cc
  Logger.cc
  Logger.h
  detail/easylogging++.h
  tmpfile_t.h
  tmpfile_t.cc
  std_ext.h
  Format.h
  detail/format.cc
  detail/format.h
)

add_library(base ${SRCS})
