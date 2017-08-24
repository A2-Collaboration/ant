#include "detail/TPCSim_tools.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"

#include "TSystem.h"
#include "TRint.h"

#include <iostream>

using namespace std;
using namespace ant;

static volatile bool interrupt = false;

int main(const int argc, const char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true;} );

    TCLAP::CmdLine cmd("TPCSim", ' ', "1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }



    return EXIT_SUCCESS;
}
