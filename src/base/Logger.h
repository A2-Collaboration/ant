#pragma once

// STL container printout made easy...
#define ELPP_STL_LOGGING
#define ELPP_DISABLE_DEFAULT_CRASH_HANDLING
#define ELPP_NO_DEFAULT_LOG_FILE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif

#include "detail/easylogging++.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#pragma GCC diagnostic pop

void SetupLogger(int argc, char* argv[]);
void SetupLogger();


// Please see
// https://github.com/easylogging/easyloggingpp
// what this nice logger can do

struct dbg {
    static int buffer_n;
};
