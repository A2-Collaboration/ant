include_directories(.)

set(BASE_ALL
  types.cc
  types.h
  interval.h
  printable.h
  printable.cc
  easylogging++.h
  Logger.cc
  Logger.h
)

add_library(base ${BASE_ALL})
