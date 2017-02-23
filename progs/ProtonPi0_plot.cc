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

    TH1D* h_IM_2g;
    TH1D* h_nPhotonsCB;

    Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t) {
        BinSettings bins_IM(100,50,250);
        h_IM_2g = HistFac.makeTH1D("IM_2g","IM / MeV","",bins_IM,"h_IM_2g");
        h_nPhotonsCB = HistFac.makeTH1D("nPhotonsCB","","",BinSettings(3),"h_nPhotonsCB");
    }

    void Fill(const Fill_t& f) const {
        h_IM_2g->Fill(f.Tree.IM_2g, f.Weight());
        h_nPhotonsCB->Fill(f.Tree.nPhotonsCB, f.Weight());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"nPhotonsCB==2", [] (const Fill_t& f) { return f.Tree.nPhotonsCB==2; } },
                              {"-", [] (const Fill_t&) { return true; } },

                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"ProtonTheta>40", [] (const Fill_t& f) { return f.Tree.Proton_Theta>40; } },
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

        cuttree::Fill<Hist_t>(myCuttree, {tree});

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
