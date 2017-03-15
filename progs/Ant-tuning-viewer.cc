#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRint.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/detail/TH2Storage.h"
#include <iostream>
#include <cstring>
#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>(
                "v",
                "verbose",
                "Verbosity level (0..9)",
                false, 0, "level");

    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>(
                    "b",
                    "batch",
                    "Run in batch mode (no GUI, autosave)",
                    false);

    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>(
                 "s",
                 "setup",
                 "Setup name",
                 false,
                 "",
                 "setup");

    auto cmd_calibration = cmd.add<TCLAP::ValueArg<string>>(  "c", "calibration","Calibration String",         true, "", "calibration");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());


    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-mcsmearing",&fake_argc,fake_argv,nullptr,0,true);


    const auto setup_name = cmd_setupname->getValue();
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    auto manager = setup->GetCalibrationDataManager();
    manager->SetOverrideToDefault(true);

    const auto id = TID(0,0,{TID::Flags_t::MC});

    const string calName = cmd_calibration->getValue();

    TCalibrationData data;
    const auto data_found = manager->GetData(calName, id, data);

    if(!data_found) {
        LOG(ERROR) << "No cdata found for " << calName;
        return EXIT_FAILURE;
    }

    TH2* h = calibration::detail::TH2Storage::Decode(data);

    canvas(calName) << drawoption("colz") << h << endc;


    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();
    setup = nullptr;
    manager = nullptr;


    return EXIT_SUCCESS;
}
