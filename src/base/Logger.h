#ifndef LOGGER_H
#define LOGGER_H

// STL container printout made easy...
#define ELPP_STL_LOGGING

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include "detail/easylogging++.h"
#pragma GCC diagnostic pop


void SetupLogger(int argc, char* argv[]);
void SetupLogger();


// Please see
// https://github.com/easylogging/easyloggingpp
// what this nice logger can do

#endif // LOGGER_H
