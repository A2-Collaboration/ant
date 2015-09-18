#include "calibration/gui/EditorWindow.h"

#include "base/Logger.h"
#include "base/CmdLine.h"

#include "TRint.h"
#include "TROOT.h"

#include <iostream>
#include <cstring>

using namespace std;
using namespace ant;
using namespace  ant::calibration::gui;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib-editor", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");

    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledValueArg<string>>("directory","calibration directory",true,"","directory");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());



    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv);

    new EditorWindow( cmd_inputfiles->getValue() );
    app->Run(kTRUE);
}
