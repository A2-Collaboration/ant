#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/etaprime/etaprime_ept.h"

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
    using Fill_t = physics::EtapEPT::Tree_t;

    TH2D* h_IM_2g;
    TH2D* h_IM_2g_wide;
    std::vector<TH2D*> h_IM_2g_Ch;

    TH2D* h_BeamPull;
    TH2D* h_BeamPull_wide;
    std::vector<TH2D*> h_BeamPull_Ch;

    TH1D* h_KinFitProb;

    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t)
    {
        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();

        HistogramFactory HistFac_IM_2g("IM_2g_Ch", HistFac);
        HistogramFactory HistFac_BeamPull("BeamPull_Ch", HistFac);

        BinSettings bins_IM(500,0,1100);
        BinSettings bins_pull(100,-5,5);

        for(unsigned ch=0;ch<ept->GetNChannels();ch++) {
            h_IM_2g_Ch.push_back(
                        HistFac_IM_2g.makeTH2D(
                            std_ext::formatter() << "IM 2g Ch=" << ch,
                            "IM / MeV","",
                            bins_IM,
                            BinSettings(ept->GetNChannels()),
                            "h_IM_2g_Ch"+to_string(ch))
                        );
            h_BeamPull_Ch.push_back(
                        HistFac_BeamPull.makeTH2D(
                            std_ext::formatter() << "Beam Pull Ch=" << ch,
                            "Pull","",
                            bins_pull,
                            BinSettings(ept->GetNChannels()),
                            "h_BeamPull_Ch"+to_string(ch))
                        );
        }

        h_IM_2g = HistFac.makeTH2D("IM 2g",
                                   "IM / MeV","",
                                   bins_IM,
                                   BinSettings(ept->GetNChannels()),
                                   "h_IM_2g");
        h_IM_2g_wide = HistFac.makeTH2D("IM 2g",
                                   "IM / MeV","",
                                   bins_IM,
                                   BinSettings(ept->GetNChannels()),
                                   "h_IM_2g_wide");

        h_BeamPull = HistFac.makeTH2D("IM 2g",
                                   "IM / MeV","",
                                   bins_pull,
                                   BinSettings(ept->GetNChannels()),
                                   "h_BeamPull");
        h_BeamPull_wide = HistFac.makeTH2D("IM 2g",
                                   "IM / MeV","",
                                   bins_pull,
                                   BinSettings(ept->GetNChannels()),
                                   "h_BeamPull_wide");

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",BinSettings(200,0,1),"h_KinFitProb");

    }


    void Fill(const Fill_t& f) const
    {
        h_IM_2g->Fill(f.IM_2g, f.TaggCh, f.TaggW);
        h_IM_2g_wide->Fill(f.IM_2g, f.TaggCh, f.TaggW_wide);
        h_IM_2g_Ch[f.TaggCh]->Fill(f.IM_2g, f.TaggCh_, f.TaggW);

        h_BeamPull->Fill(f.KinFitBeamEPull, f.TaggCh, f.TaggW);
        h_BeamPull_wide->Fill(f.KinFitBeamEPull, f.TaggCh, f.TaggW_wide);
        h_BeamPull_Ch[f.TaggCh]->Fill(f.KinFitBeamEPull, f.TaggCh_, f.TaggW);

        h_KinFitProb->Fill(f.KinFitProb, f.TaggW);
    }

    // Sig and Ref channel (can) share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"DiscardedEk=0", [] (const Fill_t& f) { return f.DiscardedEk == 0; } },
                              {"DiscardedEk<50", [] (const Fill_t& f) { return f.DiscardedEk < 50; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb>0.01", [] (const Fill_t& f) { return f.KinFitProb>0.01; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"IM Pi0", [] (const Fill_t& f) { return ParticleTypeDatabase::Pi0.GetWindow(30).Contains(f.IM_2g); } },
                              {"IM Eta", [] (const Fill_t& f) { return ParticleTypeDatabase::Eta.GetWindow(50).Contains(f.IM_2g); } },
                              {"IM EtaP", [] (const Fill_t& f) { return ParticleTypeDatabase::EtaPrime.GetWindow(50).Contains(f.IM_2g); } },
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

    physics::EtapEPT::Tree_t tree;
    {
        TTree* t;
        if(!input.GetObject("EtapEPT/tree",t)) {
            LOG(ERROR) << "Cannot find EtapEPT/tree in input";
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

    HistogramFactory HistFac("EtapEPT");

    auto cuttree = cuttree::Make<Hist_t>(HistFac, "EtapEPT", Hist_t::GetCuts());


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
