#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "TH2.h"
#include "root-addons/analysis_codes/TwoPi0_MCSmearing_tools.h"
#include "TROOT.h"
#include "TRint.h"

#include <iostream>
#include <cstring>

using namespace std;
using namespace ant;
using namespace  ant::calibration::gui;

TH2* GetHist(WrapTFileInput& file) {
    TH2* h = nullptr;
    file.GetObject("ETheta/sigma", h);
    if(!h) {
        LOG(ERROR) << "ETheta/sigma not found in " << file.FileNames();
        exit(EXIT_FAILURE);
    }
    return h;
}

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no GUI, autosave)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", true, "", "setup");
    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_data = cmd.add<TCLAP::ValueArg<string>>("d","data","Histogram file for data", true, "", "data");
    auto cmd_mc = cmd.add<TCLAP::ValueArg<string>>("mc","monte-carlo","Current iterarion mc histograms", true, "", "mc");
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

    WrapTFileInput datafile(cmd_data->getValue());
    WrapTFileInput mcfile(cmd_mc->getValue());

    const auto data_width = GetHist(datafile);
    const auto mc_width   = GetHist(mcfile);



    const auto setup_name = cmd_setupname->getValue();
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    auto manager = setup->GetCalibrationDataManager();
    manager->SetOverrideToDefault(true);

    const auto id = TID(0,0,{TID::Flags_t::MC});
    TID next;
    const auto prev = dynamic_cast<TH2*>(manager->GetTObject("CB_ClusterSmearing", "energy_smearing", id, next));

    if(!prev)  {
        LOG(ERROR) << "No prev";
        exit(1);
        //TODO: do initial step
    }

    WrapTFileOutput outfile("step.root", WrapTFileOutput::mode_t::recreate, true);
    TwoPi0_MCSmearing_Tool::CalculateUpdatedSmearing(data_width,mc_width,prev);

    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();
    setup = nullptr;

}
