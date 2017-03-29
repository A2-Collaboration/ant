#include "calibration/gui/Manager.h"
#include "calibration/gui/ManagerWindow.h"
#include "calibration/gui/AvgBuffer.h"
#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/GitInfo.h"
#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"

#include "TROOT.h"
#include "TRint.h"

#include <iostream>
#include <cstring>

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib - Fit histograms and calculate new calibration parameters", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>("c","calibration","Calibration GUI module name", true, "","calibration");
    auto cmd_sgpol = cmd.add<TCLAP::ValueArg<unsigned>>("","polyorder","Polynom order for Savitzky-Golay filter (zero is moving average)", false, 4, "polorder");
    auto cmd_average = cmd.add<TCLAP::ValueArg<unsigned>>("a","average","Average length for Savitzky-Golay filter", false, 0, "length");
    auto cmd_gotoslice = cmd.add<TCLAP::ValueArg<unsigned>>("","gotoslice","Directly skip to specified slice", false, 0, "slice");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no GUI, autosave)",false);
    auto cmd_default = cmd.add<TCLAP::SwitchArg>("","default","Put created TCalibrationData to default range",false);
    auto cmd_confirmHeaderMismatch = cmd.add<TCLAP::SwitchArg>("","confirmHeaderMismatch","Confirm mismatch in Git infos in file headers and use files anyway",false);
    auto cmd_force = cmd.add<TCLAP::SwitchArg>("","force","Ignore some safety checks (you've been warned)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "", "setup");
    auto cmd_ModuleOptions = cmd.add<TCLAP::MultiArg<string>>("O","options","Options for Calibration GUI Module, key=value",false,"");

    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_default->getValue() && cmd_average->isSet()) {
        LOG(ERROR) << "Using --default and --average (leading to ranged values) makes no sense";
        return EXIT_FAILURE;
    }

    if(cmd_gotoslice->isSet() && cmd_default->isSet()) {
        LOG(ERROR) << "Using --gotoslice with --default (no slice averaging) makes no sense";
        return EXIT_FAILURE;
    }


    auto moduleOptions = make_shared<OptionsList>();
    if(cmd_ModuleOptions->isSet()) {
        for(const auto& opt : cmd_ModuleOptions->getValue()) {
            moduleOptions->SetOption(opt);
        }
    }

    for(const auto& inputfile : cmd_inputfiles->getValue()) {
        if(std_ext::string_starts_with(inputfile, "-")) {
            LOG(ERROR) << "Found '" << inputfile << "' with starting - parsed as inputfile, might be wrongly spelled option. "
                       << "Prepend ./ to use it as inputfile.";
            return EXIT_FAILURE;
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


    unique_ptr<calibration::gui::AvgBuffer_traits<TH1>> buffer;
    if(cmd_default->isSet()) {
        buffer = std_ext::make_unique<AvgBuffer_Sum<TH1>>();
    }
    else if(cmd_average->isSet()) {
        buffer = std_ext::make_unique<AvgBuffer_SavitzkyGolay<TH1>>(
                     cmd_average->getValue(), cmd_sgpol->getValue()
                     );
    }
    if(!buffer) {
        LOG(ERROR) << "Please specify either --default or --average";
        return EXIT_FAILURE;
    }

    Manager manager(
                cmd_inputfiles->getValue(),
                move(buffer),
                cmd_confirmHeaderMismatch->getValue()
                );

    // try to find the requested calibration modules
    // the gui manager already scanned the files and provides a hint
    // for the SetupName

    const auto setup_name = cmd_setupname->isSet() ? cmd_setupname->getValue() : manager.SetupName;
    ExpConfig::Setup::SetByName(setup_name);
    auto& setup = ExpConfig::Setup::Get();

    if(cmd_default->isSet()) {
        setup.GetCalibrationDataManager()->SetOverrideToDefault(true);
    }

    GitInfo gitinfo_db(setup.GetCalibrationDataManager()->GetCalibrationDataFolder());
    if(!cmd_force->isSet() && gitinfo_db.IsDirty()) {
        LOG(ERROR) << "Cannot write to dirty calibration database without --force";
        return EXIT_FAILURE;
    }

    string calibrationguiname = cmd_calibration->getValue();
    stringstream ss_calibrationguis;
    unique_ptr<CalibModule_traits> calibrationgui = nullptr;
    for(const auto& calibration : setup.GetCalibrations()) {
        list< unique_ptr<CalibModule_traits> > guimodules;
        calibration->GetGUIs(guimodules, moduleOptions);
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
        LOG(INFO) << "Available calibrations GUIs (for setup '" << manager.SetupName << "'): "
                  << ss_calibrationguis.str();
        LOG(ERROR) << "No calibration GUI module found for given name '"
                   << calibrationguiname << "'";
        return EXIT_FAILURE;
    }

    manager.SetModule(move(calibrationgui));

    int gotoslice = cmd_gotoslice->isSet() ? cmd_gotoslice->getValue() : -1;

    if(!manager.DoInit(gotoslice)) {
        LOG(ERROR) << "Cannot initialize the calibration. Check previous messages.";
        return EXIT_FAILURE;
    }

    if(cmd_batchmode->isSet()) {
        gROOT->SetBatch();
    }

    new ManagerWindow(manager);
    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();

    return EXIT_SUCCESS;
}
