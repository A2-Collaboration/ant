include_directories(.)

set(BASE_ALL
  types.cc
  types.h
  interval.h
  printable.h
  printable.cc
  Logger.cc
  Logger.h
  detail/easylogging++.h
)

add_library(base ${BASE_ALL})
