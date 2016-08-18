#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/etaprime/etaprime_sergey.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/ProgressCounter.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile bool interrupt = false;

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct Hist_t {
    using Fill_t = physics::EtapSergey::Tree_t;


    struct IM_2g_t {

        IM_2g_t(HistogramFactory HistFac,
                BinSettings bins_IM) {

            auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();

            BinSettings bins_ThetaProton(40, 0, 40);
            BinSettings bins_ThetaPhotons(80, 0, 120);

            h = HistFac.makeTH1D(
                    "IM 2g",
                    "IM / MeV","",
                    bins_IM,
                    "h");

            h_prompt = HistFac.makeTH1D(
                           "IM 2g",
                           "IM / MeV","",
                           bins_IM,
                           "h_prompt");
            h_random = HistFac.makeTH1D(
                           "IM 2g",
                           "IM / MeV","",
                           bins_IM,
                           "h_random");

            h_TaggCh = HistFac.makeTH2D(
                           "IM 2g",
                           "IM / MeV","TaggCh",
                           bins_IM,
                           BinSettings(ept->GetNChannels()),
                           "h_TaggCh");

            h_ThetaProton = HistFac.makeTH2D(
                                "IM 2g",
                                "IM / MeV","#theta / #circ",
                                bins_IM,
                                bins_ThetaProton,
                                "h_ThetaProton");

            h_ThetaPhotons = HistFac.makeTH2D(
                                 "IM 2g",
                                 "IM / MeV","#theta / #circ",
                                 bins_IM,
                                 bins_ThetaPhotons,
                                 "h_ThetaPhotons");
        }

        void Fill(double IM_2g, const Fill_t& f) const {
            h->Fill(IM_2g, f.TaggW);
            h_prompt->Fill(IM_2g, f.TaggW>0);
            h_random->Fill(IM_2g, -(f.TaggW<0)*f.TaggW);
            h_TaggCh->Fill(IM_2g, f.TaggCh, f.TaggW);
            h_ThetaProton->Fill(IM_2g, f.ProtonTheta, f.TaggW);
            for(const auto& PhotonTheta : f.PhotonsTheta()) {
                h_ThetaPhotons->Fill(IM_2g, PhotonTheta, f.TaggW);
            }
        }

        TH1D* h;
        TH1D* h_prompt;
        TH1D* h_random;
        TH2D* h_TaggCh;
        TH2D* h_ThetaProton;
        TH2D* h_ThetaPhotons;
    };


    const BinSettings bins_IM_etap{150, 800, 1050};
    const BinSettings bins_IM_eta {150, 400,  700};

    TH1D* h_KinFitProb;
    TH2D* h_MissingMass;

    IM_2g_t IM_2g;
    IM_2g_t IM_2g_raw;
    IM_2g_t IM_2g_eta;
    IM_2g_t IM_2g_raw_eta;

    TH2D* h_IM_2g_raw_fitted;

    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t) :
        IM_2g(HistogramFactory("IM_2g",HistFac), bins_IM_etap),
        IM_2g_raw(HistogramFactory("IM_2g_raw",HistFac), bins_IM_etap),
        IM_2g_eta(HistogramFactory("IM_2g_eta",HistFac), bins_IM_eta),
        IM_2g_raw_eta(HistogramFactory("IM_2g_raw_eta",HistFac), bins_IM_eta)
    {
        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();


        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",BinSettings(200,0,1),"h_KinFitProb");

        BinSettings bins_MM(150,790,1090);

        h_MissingMass = HistFac.makeTH2D("Missing Mass",
                                         "MM / MeV","",
                                         bins_MM,
                                         BinSettings(ept->GetNChannels()),
                                         "h_MissingMass");

        h_IM_2g_raw_fitted = HistFac.makeTH2D("IM 2#gamma raw vs. fitted","IM raw","IM fitted",
                                              bins_IM_etap, bins_IM_etap, "h_IM_2g_raw_fitted");
    }


    void Fill(const Fill_t& f) const
    {
        IM_2g.Fill(f.FittedPhotonSum, f);
        IM_2g_raw.Fill(f.PhotonSum, f);
        IM_2g_eta.Fill(f.FittedPhotonSum, f);
        IM_2g_raw_eta.Fill(f.PhotonSum, f);

        h_IM_2g_raw_fitted->Fill(f.PhotonSum, f.FittedPhotonSum, f.TaggW);

        h_KinFitProb->Fill(f.KinFitProb, f.TaggW);
        h_MissingMass->Fill(f.MissingMass, f.TaggCh, f.TaggW);
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"+3<ZVertex<+5", [] (const Fill_t& f) { return f.TrueZVertex>+3 && f.TrueZVertex<+5; } },
                              {"-1<ZVertex<+1", [] (const Fill_t& f) { return f.TrueZVertex>-1 && f.TrueZVertex<+1; } },
                              {"-5<ZVertex<-3", [] (const Fill_t& f) { return f.TrueZVertex>-5 && f.TrueZVertex<-3; } },
                          });

        return cuts;
    }
};

int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true; } );

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "Setup_2014_07_EPT_Prod", "setup");

    cmd.parse(argc, argv);

    const auto setup_name = cmd_setupname->getValue();
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }



    WrapTFileInput input(cmd_input->getValue());

    physics::EtapSergey::Tree_t tree;
    {
        TTree* t;
        if(!input.GetObject("EtapSergey/tree",t)) {
            LOG(ERROR) << "Cannot find EtapSergey/tree in input";
            return 1;
        }
        if(!tree.Matches(t)) {
            LOG(ERROR) << "Found tree does not match";
            return 1;
        }
        tree.LinkBranches(t);
    }


    auto entries = tree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    HistogramFactory HistFac("EtapSergey");

    auto cuttree = cuttree::Make<Hist_t>(HistFac, "EtapSergey", Hist_t::GetCuts());


    LOG(INFO) << "Tree entries=" << entries;
    auto max_entries = entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    long long entry = 0;
    ProgressCounter::Interval = 3;
    ProgressCounter progress(
                [&entry, entries] (std::chrono::duration<double>) {
        LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    });

    for(entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        tree.Tree->GetEntry(entry);

        cuttree::Fill<Hist_t>(cuttree, tree);

        ProgressCounter::Tick();
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("EtapOmegaG_plot",&argc,argv,nullptr,0,true);

            if(masterFile)
                LOG(INFO) << "Stopped running, but close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }


    return 0;
}
