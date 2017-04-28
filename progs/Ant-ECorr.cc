#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "tclap/CmdLine.h"
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
#include "base/TH_ext.h"
#include <array>
#include "analysis/plot/RootDraw.h"
#include "base/std_ext/math.h"
#include "analysis/plot/HistogramFactory.h"
#include "base/Array2D.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::calibration::gui;

TH2* GetHist(WrapTFileInput& file, const string& histname) {
    TH2* h = nullptr;
    file.GetObject(histname, h);
    if(!h) {
        LOG(ERROR) << histname << " not found in " << file.FileNames();
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

struct TCLAPBinSettings : BinSettings {
    using BinSettings::BinSettings;
    using ValueCategory = TCLAP::ValueLike;
};

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>   ("v","verbose", "Verbosity level (0..9)", false, 0,"level");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>       ("b","batch",   "Run in batch mode (no GUI, autosave)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("", "write-to-setup",   "Store ECorr of this setup in the database. Do not store if not set.",       false, "", "setup");
    auto cmd_detector  = cmd.add<TCLAP::ValueArg<string>>("" ,"detector","Detector Name",    true, "", "detector");
    auto cmd_file      = cmd.add<TCLAP::ValueArg<string>>("" ,"file",    "Input file",       true, "", "file");
    auto cmd_bins      = cmd.add<TCLAP::ValueArg<TCLAPBinSettings>>("" ,"bins",    "Input file",       false, TCLAPBinSettings(0), "file");
    auto cmd_userecovery = cmd.add<TCLAP::SwitchArg>       ("", "use-recovery",   "Plot energy recovery instead of correction factor. cannot be used together with write-to-setup", false);

    cmd.parse(argc, argv);
    const bool SaveToDatabase = cmd_setupname->isSet();

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    const auto det = ToUpper(cmd_detector->getValue());

    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-mcsmearing",&fake_argc,fake_argv,nullptr,0,true);


    const string histname = std_ext::formatter() << "MCClusterECorr/" << det << "/h_EtrueErec";
    const string statsHistName = std_ext::formatter() << "MCClusterECorr/" << det << "/h_nFills";

    const auto invert = [&cmd_userecovery] (TH2* h) {
        if(cmd_userecovery->isSet()) {
            for(int x=1;x<=h->GetNbinsX();++x) {
                for(int y=1;y<=h->GetNbinsY();++y) {
                    h->SetBinContent(x,y,1./h->GetBinContent(x,y));
                }
            }
        }
        return h;
    };

    WrapTFileInput infile(cmd_file->getValue());
    const auto h_ecorr  = invert(GetHist(infile, histname));
    const auto h_nfills = GetHist(infile, statsHistName);

    const auto range = interval<double>::CenterWidth(1.0,.15);

    std_ext::IQR iqr;

    for(int x=1;x<=h_ecorr->GetNbinsX();++x) {
        for(int y=1;y<=h_ecorr->GetNbinsY();++y) {
             if(h_nfills->GetBinContent(x,y) > 500.0) {
                 const auto v = h_ecorr->GetBinContent(x,y);
                 iqr.Add(v);
             }
        }
    }

    analysis::HistogramFactory f("ECorr");

    auto h       = f.makeTH1D("Factors",           "ECorr Factor","", BinSettings(250,0.8,1.05));

    for(int x=1;x<=h_ecorr->GetNbinsX();++x) {
        for(int y=1;y<=h_ecorr->GetNbinsY();++y) {
             if(h_nfills->GetBinContent(x,y) > 500.0) {
                 const auto v = h_ecorr->GetBinContent(x,y);
                     h->Fill(v);
             }
        }
    }

    const auto TitleAppend = [] (TH2* h, const string& s) {
        const string t = h->GetTitle() + s;
        h->SetTitle(t.c_str());
    };

    const auto norm = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) - 1.0;
    LOG(INFO) << norm;

    auto stat_cutted = static_cast<TH2D*>(TH_ext::Apply(h_ecorr, h_nfills, [range, norm] (const double ecorr, const double n) {
        return n > 500.0 ? ecorr : std_ext::NaN;
    }));

    TitleAppend(stat_cutted, ": n>500");


    auto shifted = static_cast<TH2D*>(TH_ext::Apply(stat_cutted, [range, norm] (const double ecorr) {
        return range.Clip(ecorr-norm);
    }));
    TitleAppend(shifted, ": shifted");

    auto filled = TH_ext::Clone(shifted, "ECorrFilled");
    TitleAppend(filled, ": filled");

    Array2D_TH2D a(filled);

    a.FloodFillAverages();


    canvas(formatter() << "ECorr: " << cmd_file->getValue())
            << drawoption("colz")
            << h << stat_cutted << shifted << filled
            << endc;


    if(SaveToDatabase) {
        LOG(INFO) << "Writing to calibration database";
        const auto setup_name = cmd_setupname->getValue();
        ExpConfig::Setup::SetByName(setup_name);
        auto& setup = ExpConfig::Setup::Get();

        auto manager = setup.GetCalibrationDataManager();
        manager->SetOverrideToDefault(true);

        const auto id = TID(0,0,{TID::Flags_t::MC});


        const string calName = std_ext::formatter() << det << "_ClusterECorr";

        TCalibrationData cdata(calName, id, id);
        calibration::detail::TH2Storage::Encode(filled, cdata);

        manager->Add(cdata,  Calibration::AddMode_t::AsDefault);
    }


    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();

    return EXIT_SUCCESS;
}
