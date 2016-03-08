#pragma once

// STL container printout made easy...
#define ELPP_STL_LOGGING
#define ELPP_DISABLE_DEFAULT_CRASH_HANDLING
#define ELPP_NO_DEFAULT_LOG_FILE
#define ELPP_THREAD_SAFE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wunused-parameter"


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

namespace ant {
namespace logger {

struct DebugInfo {
    static int nUnpackedBuffers;
    static long long nProcessedEvents;
};

}}

