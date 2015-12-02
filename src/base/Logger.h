#pragma once

// STL container printout made easy...
#define ELPP_STL_LOGGING
#define ELPP_DISABLE_DEFAULT_CRASH_HANDLING
#define ELPP_NO_DEFAULT_LOG_FILE

#pragma GCC diagnostic push
#pragma clang diagnostic push

#pragma GCC diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wpessimizing-move"

#include "detail/easylogging++.h"

#pragma GCC diagnostic pop
#pragma clang diagnostic pop

void SetupLogger(int argc, char* argv[]);
void SetupLogger();


// Please see
// https://github.com/easylogging/easyloggingpp
// what this nice logger can do
