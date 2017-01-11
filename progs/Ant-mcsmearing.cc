#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "TH2.h"
#include "TH3.h"
#include "root-addons/analysis_codes/TwoPi0_MCSmearing_tools.h"
#include "TROOT.h"
#include "TRint.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/detail/TH2Storage.h"
#include "base/TH_ext.h"
#include <iostream>
#include <cstring>
#include "TDirectory.h"
#include "base/std_ext/string.h"

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

inline string ToUpper(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::toupper);
    return out;
}

inline string ToLower(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::tolower);
    return out;
}

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>("b","batch","Run in batch mode (no GUI, autosave)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup name", false, "", "setup");

    auto cmd_detector  = cmd.add<TCLAP::ValueArg<string>>("","detector","Detector Name",         true, "", "detector");

    auto cmd_fit  = cmd.add<TCLAP::ValueArg<string>>("f","fit","Fit peaks in IM spectra",         false, "", "fit");
    auto cmd_step = cmd.add<TCLAP::SwitchArg>("","step", "Calculate next iteration step. Init if no previous data in calib database", false);
    auto cmd_data = cmd.add<TCLAP::ValueArg<string>>("d","data","Histogram file for data",         false, "", "data");
    auto cmd_mc   = cmd.add<TCLAP::ValueArg<string>>("m","mc",  "Current iterarion mc histograms", false, "", "mc");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(!cmd_detector->isSet()) {
        LOG(ERROR) << "Detector not set";
        return EXIT_FAILURE;
    }

    const auto det = ToLower(cmd_detector->getValue());

    cout << det << " " << cmd_detector->getValue() << endl;

    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-mcsmearing",&fake_argc,fake_argv,nullptr,0,true);

    if(cmd_fit->isSet() && !cmd_step->isSet()) {
        // FIT
        TH3* h_ETheta = nullptr;
        WrapTFileInput infile(cmd_fit->getValue());
        const string histname = std_ext::formatter() << "TwoPi0_MCSmearing/m2Pi0/" << det << "_pi0_ETheta";
        infile.GetObject(histname, h_ETheta);

        if(!h_ETheta) {
            LOG(ERROR) << "Histogram " << histname << " not found!" << endl;
            return EXIT_FAILURE;
        }

        WrapTFileOutput f( std_ext::formatter() << cmd_fit->getValue() << "." << det << "_fitted.root", WrapTFileOutput::mode_t::recreate, true);

        auto d = gDirectory->mkdir("ETheta");
        d->cd();
        ant::TwoPi0_MCSmearing_Tool::AnalyseChannelE(h_ETheta);

        app->Run(kTRUE);
    }
    else if(!cmd_fit->isSet() && cmd_step->isSet()) {

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
        const string calName = std_ext::formatter() << ToUpper(det) << "_ClusterSmearing";
        TID next;
        TCalibrationData prev_data;
        const auto prev_avail = manager->GetData(calName, id, prev_data, next);

        TH2* smearing = nullptr;

        if(prev_avail) {
            TH2* prev_hist = calibration::detail::TH2Storage::Decode(prev_data);
            smearing = TwoPi0_MCSmearing_Tool::CalculateUpdatedSmearing(data_width,mc_width, prev_hist);
        } else {
            smearing = TwoPi0_MCSmearing_Tool::CalculateInitialSmearing(data_width, mc_width);
        }

        TCalibrationData cdata(calName, id, id);
        calibration::detail::TH2Storage::Encode(smearing, cdata);

        manager->Add(cdata,  Calibration::AddMode_t::AsDefault);

        app->Run(kTRUE);
        ExpConfig::Setup::Cleanup();
        setup = nullptr;
        manager = nullptr;
    }
    else {
        LOG(ERROR) << "Don't specify --fit and --step at the same time. First fit, then step.";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
