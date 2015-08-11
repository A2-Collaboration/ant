
#include "base/std_ext.h"

#include "calibration/gui/debug_GUIModule.h"
#include "calibration/gui/Manager.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"

#include "TRint.h"
#include "TExec.h"

#include <iostream>

using namespace std;
using namespace ant;

using namespace  ant::calibration::gui;

struct MyExec : TExec {

    const TRint* rint = nullptr;
    Manager* gui;

    MyExec(Manager* gui_, const TRint* Rint) :
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


    TCLAP::CmdLine cmd("ant_gui", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration GUI module name", true, "","modulename");
    auto cmd_averagelength = cmd.add<TCLAP::ValueArg<int>>("a","average","Average length for moving window (zero sums everything up)", false, 0, "length");
    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    unique_ptr<Manager> gui = std_ext::make_unique<Manager>(
                                  cmd_inputfiles->getValue(),
                                  cmd_averagelength->getValue()
                                  );

    // try to find the requested calibration modules
    // the gui manager already scanned the files and provides a hint
    // for the SetupName

    auto setup = ExpConfig::Setup::Get(gui->SetupName);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name '" << gui->SetupName << "' (extracted from given files)";
        return 1;
    }

    string calibrationguiname = cmd_calibration->getValue();
    stringstream ss_calibrationguis;
    shared_ptr<Manager_traits> calibrationgui = nullptr;
    for(const auto& calibration : setup->GetCalibrations()) {
        list< unique_ptr<Manager_traits> > guimodules;
        calibration->GetGUIs(guimodules);
        for(auto& guimodule : guimodules) {
            ss_calibrationguis << guimodule->GetName() << " ";
            if(guimodule->GetName() == calibrationguiname) {
                if(calibrationgui != nullptr) {
                    LOG(ERROR) << "Found two calibration GUI modules with name '"
                               << calibrationguiname << "'";
                    return 1;
                }
                calibrationgui = move(guimodule);
            }
        }
    }
    if(calibrationgui == nullptr) {
        LOG(INFO) << "Available calibrations GUIs (for setup '" << gui->SetupName << "''): "
                  << ss_calibrationguis.str();
        LOG(ERROR) << "No calibration GUI module found for given name '"
                   << calibrationguiname << "'";
        return 1;
    }

    gui->SetModule(calibrationgui);

    int fake_argc=0;
    auto app = std_ext::make_unique<TRint>("ant_gui",&fake_argc,nullptr);

    auto exec = std_ext::make_unique<MyExec>(gui.get(), app.get());
    exec->Exec("firstcall");
    app->Run(kTRUE);

    ExpConfig::Setup::Cleanup();
    setup = nullptr;
    gui = nullptr;
    calibrationgui = nullptr;
    app = nullptr;
    exec = nullptr;
}
