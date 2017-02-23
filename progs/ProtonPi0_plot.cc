#include "analysis/physics/pi0/ProtonPi0.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"
#include "expconfig/ExpConfig.h"
#include "base/ProgressCounter.h"
#include "analysis/plot/CutTree.h"

#include "TRint.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;


struct Hist_t {
    using Tree_t = physics::ProtonPi0::Tree_t;
    struct Fill_t {
        Fill_t(const Tree_t& t) : Tree(t) {}
        const Tree_t& Tree;
        double Weight() const {
            return Tree.TaggW;
        }
    };

    TH1D* h_KinFitProb;
    TH1D* h_TaggT;
    TH1D* h_IM_2g;
    TH1D* h_IM_2g_fitted;
    TH1D* h_nPhotonsCB;
    TH2D* h_ProtonThetaEk;
    TH2D* h_ProtonMinPIDPhi;
    TH1D* h_ProtonMinPIDCh;
    TH1D* h_ProtonPhi;
    TH2D* h_Banana;


    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t) {
        BinSettings bins_IM(100,50,250);

        auto pid = ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID);

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",BinSettings(100,0,1),"h_KinFitProb");
        h_TaggT = HistFac.makeTH1D("TaggT-CBAvgTime","t / ns","",BinSettings(200,-40,40),"h_TaggT");
        h_IM_2g = HistFac.makeTH1D("IM_2g","IM / MeV","",bins_IM,"h_IM_2g");
        h_IM_2g_fitted = HistFac.makeTH1D("IM_2g_fitted","IM / MeV","",bins_IM,"h_IM_2g_fitted");
        h_nPhotonsCB = HistFac.makeTH1D("nPhotonsCB","","",BinSettings(3),"h_nPhotonsCB");
        h_ProtonThetaEk = HistFac.makeTH2D("Proton #theta E_k","Ek / MeV","#theta / #circ",
                                           BinSettings(100,0,1000),BinSettings(40,20,80),
                                           "h_ProtonThetaEk");
        h_ProtonMinPIDPhi = HistFac.makeTH2D("Proton MinPID","#Delta#phi / #circ","PID Channel",
                                             BinSettings(200,-200,200),BinSettings(pid->GetNChannels()),
                                             "h_ProtonMinPIDPhi");
        h_ProtonMinPIDCh = HistFac.makeTH1D("Proton MinPIDCh","PID Channel","",
                                            BinSettings(pid->GetNChannels()),
                                            "h_ProtonMinPIDCh");
        h_ProtonPhi = HistFac.makeTH1D("Proton Phi","#phi / #circ","",BinSettings(pid->GetNChannels(), -180, 180),"h_ProtonPhi");
        h_Banana = HistFac.makeTH2D("Banana","E_k / MeV","VetoE / MeV",
                                    BinSettings(100,0,1000),BinSettings(100,0,10),
                                    "h_Banana");
    }

    void Fill(const Fill_t& f) const {
        h_KinFitProb->Fill(f.Tree.FitProb);
        h_TaggT->Fill(f.Tree.TaggT-f.Tree.CBAvgTime);
        h_IM_2g->Fill(f.Tree.IM_2g);
        h_IM_2g_fitted->Fill(f.Tree.IM_2g_fitted);
        h_nPhotonsCB->Fill(f.Tree.nPhotonsCB);
        h_ProtonThetaEk->Fill(f.Tree.Proton_Ek, f.Tree.Proton_Theta);
        h_ProtonMinPIDPhi->Fill(f.Tree.Proton_MinPIDPhi,f.Tree.Proton_MinPIDCh);
        h_ProtonMinPIDCh->Fill(f.Tree.Proton_MinPIDCh);
        h_ProtonPhi->Fill(f.Tree.Proton_Phi);
        h_Banana->Fill(f.Tree.Proton_Ek, f.Tree.Proton_VetoE);
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        using i_t = interval<double>;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"prompt", [] (const Fill_t& f) { return i_t(-2.5, 2.5).Contains(f.Tree.TaggT-f.Tree.CBAvgTime); }},
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"nPhotonsCB>0", [] (const Fill_t& f) { return f.Tree.nPhotonsCB>0; } },
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"100<IM<180", [] (const Fill_t& f) { return i_t(100, 180).Contains(f.Tree.IM_2g_fitted); }},
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"ProtonTheta>40", [] (const Fill_t& f) { return f.Tree.Proton_Theta>40; } },
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"-25<Phi<25", [] (const Fill_t& f) { return i_t(-25,25).Contains(f.Tree.Proton_MinPIDPhi); } },
                              {"-", [] (const Fill_t&) { return true; } },
                          });
        return cuts;
    }
};

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup","Override setup name", false, "Setup_2014_07_EPT_Prod", "setup");
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    cmd.parse(argc, argv);

    const auto setup_name = cmd_setupname->getValue() ;
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    WrapTFileInput input(cmd_input->getValue());

    TTree* t = nullptr;
    if(!input.GetObject("ProtonPi0/t",t)) {
        LOG(ERROR) << "Cannot find tree in input file";
    }

    analysis::physics::ProtonPi0::Tree_t tree;
    tree.LinkBranches(t);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }

    HistogramFactory histFac("ProtonPi0");

    auto myCuttree = cuttree::Make<Hist_t>(histFac, "tree", Hist_t::GetCuts());

    auto max_entries = tree.Tree->GetEntries();

    LOG(INFO) << "Max tree entries=" << max_entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<tree.Tree->GetEntries()) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }


    for(long long entry=0;entry<max_entries;entry++) {
        tree.Tree->GetEntry(entry);

        cuttree::Fill<Hist_t>(myCuttree, tree);

    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            argc=0; // prevent TRint to parse any cmdline
            TRint app("ProtonPi0_plot",&argc,argv,nullptr,0,true);

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
