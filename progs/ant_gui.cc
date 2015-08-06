
#include "base/std_ext.h"

#include "TRint.h"
#include "TCanvas.h"
#include "TThread.h"
#include "RQ_OBJECT.h"
#include "TExec.h"
#include "TObject.h"
#include "base/Logger.h"
#include "calibration/gui/debug_GUIModule.h"
#include "calibration/gui/manager/CalibrationGUI.h"
#include "calibration/gui/CalCanvas.h"
#include "base/CmdLine.h"

#include <iostream>

using namespace std;
using namespace ant;

using namespace  ant::calibration::gui;

struct MyExec : TExec {

    const TRint* rint = nullptr;
    CalibrationGUI* gui;

    MyExec(CalibrationGUI* gui_, const TRint* Rint) :
        rint(Rint),
        gui(gui_)
    {
        gui->ConnectReturnFunc("TExec", this, "Exec(=\"\")");
    }

    virtual void Exec(const char* arg) override {
        if(string(arg) != "firstcall" && !rint->IsRunning()) {
            return;
        }
        while(gui->Run()) {}
    }
};

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("ant", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_input  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
    el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    auto d = std_ext::make_unique<ant::calibration::gui::DebugModule>();

    int fake_argc=0;
    auto app = std_ext::make_unique<TRint>("app",&fake_argc,nullptr);

    unique_ptr<CalibrationGUI> gui = std_ext::make_unique<CalibrationGUI>(move(d),5);
    gui->SetFileList(cmd_input->getValue());

    auto exec = std_ext::make_unique<MyExec>(gui.get(), app.get());
    exec->Exec("firstcall");
    app->Run(kTRUE);

    gui = nullptr;
    app = nullptr;
    exec = nullptr;
}
