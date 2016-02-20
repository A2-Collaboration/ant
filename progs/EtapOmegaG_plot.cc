#include "base/Logger.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"

#include "base/CmdLine.h"
#include "base/interval.h"

#include "base/iterators.h"
#include "base/std_ext/string.h"

#include "TH1D.h"

#include <string>
#include <list>

using namespace ant;
using namespace std;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");

    auto cmd_file = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");;

    cmd.parse(argc, argv);

    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("plot",&fake_argc,fake_argv,nullptr,0,true);


    app.Run();

    return 0;
}