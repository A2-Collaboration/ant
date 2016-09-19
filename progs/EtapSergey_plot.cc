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
#include "base/ParticleType.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile bool interrupt = false;

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)
    {
        using histstyle::Mod_t;
        this->GetHist(0, "Data", Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2.0));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 2.0));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+2));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = f.MCTrue;

        auto get_bkg_name = [] (unsigned mctrue) {
            const string& name = mctrue>=10 ?
                                     physics::EtapSergey::ptreeBackgrounds[mctrue-10].Name
                                 : "Other";
            return "Bkg_"+name;
        };

        using histstyle::Mod_t;
        const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      Mod_t::MakeLine(histstyle::color_t::Get(mctrue-9), 1, kGray+2)
                                                      );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 9)
                this->GetHist(4).Fill(f);
        }
    }
};

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct Hist_t {
    using Fill_t = physics::EtapSergey::Tree_t;

    TH1D* h_KinFitProb;
    TH2D* h_MissingMass;

    TH1D* h_IM_4g;
    TH2D* h_IM_3g_4g_high;
    TH1D* h_CBVetoSumE;
    TH2D* h_gNonPi0_CaloE_Theta_lo;
    TH2D* h_gNonPi0_CaloE_Theta_hi;


    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t)
    {
//        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();
        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",BinSettings(200,0,1),"h_KinFitProb");

        BinSettings bins_IM(100,600,1000);

        h_IM_3g_4g_high = HistFac.makeTH2D("IM 3g vs. 4g high","IM 4g","IM 3g high",bins_IM,bins_IM,"h_IM_3g_4g_high");
        h_IM_4g = HistFac.makeTH1D("IM 4g","IM 4g","",bins_IM,"h_IM_4g");

        h_gNonPi0_CaloE_Theta_lo = HistFac.makeTH2D("gNonPi0 low","Theta","CaloE",
                                                 BinSettings(180),BinSettings(100,0,300),
                                                 "h_gNonPi0_CaloE_Theta_lo");
        h_gNonPi0_CaloE_Theta_hi = HistFac.makeTH2D("gNonPi0 high","Theta","CaloE",
                                                 BinSettings(180),BinSettings(100,0,300),
                                                 "h_gNonPi0_CaloE_Theta_hi");


        h_CBVetoSumE = HistFac.makeTH1D("CB Veto SumE","VetoE / MeV","",BinSettings(100,0,20),"h_CBVetoSumE");
    }


    void Fill(const Fill_t& f) const
    {
        if(f.TaggT < -50.0 || f.TaggT > 50.0)
          return;
        const auto TaggW = f.TaggT > -5.5 && f.TaggT < 5.5 ? 1.0 : -(5.5+5.5)/(50+50-(5.5+5.5));

        h_KinFitProb->Fill(f.KinFitProb, TaggW);

        h_IM_4g->Fill(f.IM_4g, TaggW);
        h_IM_3g_4g_high->Fill(f.IM_4g, f.IM_3g()[1], TaggW);

        {
            auto& caloEs = f.gNonPi0_CaloE();
            auto i_minE = caloEs.front() < caloEs.back() ? 0 : 1;
            auto i_maxE = 1-i_minE;
            h_gNonPi0_CaloE_Theta_lo->Fill(f.gNonPi0_Theta().at(i_minE), f.gNonPi0_CaloE().at(i_minE), TaggW);
            h_gNonPi0_CaloE_Theta_hi->Fill(f.gNonPi0_Theta().at(i_maxE), f.gNonPi0_CaloE().at(i_maxE), TaggW);
        }

        h_CBVetoSumE->Fill(f.CBVetoSumE, TaggW);
    }

    std::vector<TH1*> GetHists() const {
        return {h_KinFitProb, h_IM_4g, h_CBVetoSumE};
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"Simple", [] (const Fill_t& f) { return f.TaggE>1447 && f.AntiPi0FitProb < 0.00001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"Prob", [] (const Fill_t& f) { return f.TreeFitProb>0.2 && f.AntiEtaFitProb < 0.00001; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"IM_Omega", [] (const Fill_t& f) {
                                   const auto& window = ParticleTypeDatabase::Omega.GetWindow(80);
                                   return window.Contains(f.IM_3g()[1]);
                               }},
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"CBVetoSumE", [] (const Fill_t& f) { return f.CBVetoSumE<0.4; }},
                              {"-", [] (const Fill_t& ) { return true; }},
                          });
        auto gNonPi0_cut_1 = [] (const Fill_t& f) {
            auto& caloEs = f.gNonPi0_CaloE();
            auto i_minE = caloEs.front() < caloEs.back() ? 0 : 1;
            auto theta = f.gNonPi0_Theta().at(i_minE);
            auto caloE = caloEs.at(i_minE);
            return caloE > 230.0*(1.0-theta/160.0);
        };

        auto gNonPi0_cut_2 = [] (const Fill_t& f) {
            auto& caloEs = f.gNonPi0_CaloE();
            auto i_minE = caloEs.front() < caloEs.back() ? 0 : 1;
            auto theta = f.gNonPi0_Theta().at(i_minE);
            auto caloE = caloEs.at(i_minE);
            if(theta<22)
                return caloE > 140;
            else
                return caloE > 60;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"gNonPi0_1", gNonPi0_cut_1},
                              {"gNonPi0_2", gNonPi0_cut_2},
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
        if(!input.GetObject("EtapSergey/treeSergey",t)) {
            LOG(ERROR) << "Cannot find EtapSergey/treeSergey in input";
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

    auto cuttree = cuttree::Make<MCTrue_Splitter<Hist_t>>(HistFac,
                                                          "EtapSergey",
                                                          Hist_t::GetCuts()
                                                          );


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

        cuttree::Fill<MCTrue_Splitter<Hist_t>>(cuttree, tree);

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
