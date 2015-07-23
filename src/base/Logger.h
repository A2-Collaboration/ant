#pragma once

// STL container printout made easy...
#define ELPP_STL_LOGGING
#define ELPP_DISABLE_DEFAULT_CRASH_HANDLING

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include "detail/easylogging++.h"
#pragma GCC diagnostic pop


void SetupLogger(int argc, char* argv[]);
void SetupLogger();


// Please see
// https://github.com/easylogging/easyloggingpp
// what this nice logger can do
