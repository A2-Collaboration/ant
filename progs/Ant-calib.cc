#include "calibration/gui/Manager.h"
#include "calibration/gui/ManagerWindow.h"
#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"

#include "TROOT.h"
#include "TRint.h"

#include <iostream>
#include <cstring>

using namespace std;
using namespace ant;
using namespace  ant::calibration::gui;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib - Fit histograms and calculate new calibration parameters", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration GUI module name", true, "","calibration");
    auto cmd_averagelength = cmd.add<TCLAP::ValueArg<unsigned>>("a","average","Average length for moving window (zero sums everything up)", false, 0, "length");
    auto cmd_gotoslice = cmd.add<TCLAP::ValueArg<unsigned>>("","gotoslice","Directly skip to specified slice", false, 0, "slice");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no GUI, autosave)",false);
    auto cmd_default = cmd.add<TCLAP::SwitchArg>("","default","Put created TCalibrationData to default range",false);
    auto cmd_confirmHeaderMismatch = cmd.add<TCLAP::SwitchArg>("","confirmHeaderMismatch","Confirm mismatch in Git infos in file headers and use files anyway",false);
    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_default->getValue() && cmd_averagelength->getValue()>0) {
        LOG(ERROR) << "Default and using moving window leading to ranged values makes no sense";
        return 1;
    }

    if(cmd_gotoslice->isSet() && cmd_averagelength->getValue()==0) {
        LOG(ERROR) << "Goto slice without averaging makes no sense";
        return 1;
    }

    for(auto inputfile : cmd_inputfiles->getValue()) {
        if(std_ext::string_starts_with(inputfile, "-")) {
            LOG(ERROR) << "Found '" << inputfile << "' with starting - parsed as inputfile, might be wrongly spelled option. "
                       << "Prepend ./ to use it as inputfile.";
            return 1;
        }
    }

    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv,nullptr,0,true);


    auto manager = std_ext::make_unique<Manager>(
                       cmd_inputfiles->getValue(),
                       cmd_averagelength->getValue(),
                       cmd_confirmHeaderMismatch->getValue()
                       );

    // try to find the requested calibration modules
    // the gui manager already scanned the files and provides a hint
    // for the SetupName

    auto setup = ExpConfig::Setup::Get(manager->SetupName);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name '" << manager->SetupName << "' (extracted from input files)";
        return 1;
    }

    if(cmd_default->isSet()) {
        setup->GetCalibrationDataManager()->SetOverrideToDefault(true);
    }

    string calibrationguiname = cmd_calibration->getValue();
    stringstream ss_calibrationguis;
    shared_ptr<CalibModule_traits> calibrationgui = nullptr;
    for(const auto& calibration : setup->GetCalibrations()) {
        list< unique_ptr<CalibModule_traits> > guimodules;
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
        LOG(INFO) << "Available calibrations GUIs (for setup '" << manager->SetupName << "'): "
                  << ss_calibrationguis.str();
        LOG(ERROR) << "No calibration GUI module found for given name '"
                   << calibrationguiname << "'";
        return 1;
    }

    manager->SetModule(calibrationgui);

    int gotoslice = cmd_gotoslice->isSet() ? cmd_gotoslice->getValue() : -1;

    if(!manager->DoInit(gotoslice)) {
        LOG(ERROR) << "Cannot initialize the calibration. Check previous messages.";
        return 1;
    }

    if(cmd_batchmode->isSet()) {
        gROOT->SetBatch();
    }

    new ManagerWindow(manager.get());
    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();
    setup = nullptr;
    manager = nullptr;
    calibrationgui = nullptr;
}
