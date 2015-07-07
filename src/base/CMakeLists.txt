include_directories(.)

set(BASE_ALL
    types.cc
    types.h
    interval.h
    printable.h
    printable.cc
)

add_library(base SHARED ${BASE_ALL})
list(APPEND LIBS base)
