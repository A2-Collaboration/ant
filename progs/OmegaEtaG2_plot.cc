#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/utils/particle_tools.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile bool interrupt = false;

class MyTInterruptHandler : public TSignalHandler {
public:
    MyTInterruptHandler() : TSignalHandler(kSigInterrupt, kFALSE) { }

    Bool_t  Notify() {
        if (fDelay) {
            fDelay++;
            return kTRUE;
        }
        interrupt = true;
        cout << " >>> Interrupted! " << endl;
        return kTRUE;
    }
};

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    const decltype (physics::OmegaEtaG2::makeChannels()) Channels;

    MCTrue_Splitter(const HistogramFactory& histFac) : cuttree::StackedHists_t<Hist_t>(histFac),
      Channels(physics::OmegaEtaG2::makeChannels())
    {
        using cuttree::HistMod_t;

        // TODO: derive this from channel map
        this->GetHist(0, "Data", HistMod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  HistMod_t::MakeLine(kRed, 2.0));
        this->GetHist(2, "Ref",  HistMod_t::MakeLine(kGreen, 2.0));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", HistMod_t::MakeLine(kBlack, 1.0));
        this->GetHist(4, "Bkg_MC", HistMod_t::MakeLine(kGray, 1.0));
    }

    void Fill(const Fill_t& f) {

        const int mctrue = f.Tree.Channel;

        auto get_bkg_name = [] (unsigned mctrue) {
            const auto entry = physics::OmegaEtaG2::reaction_channels.channels.find(mctrue);

            if(entry!=physics::OmegaEtaG2::reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        using cuttree::HistMod_t;
        const Hist_t& hist = mctrue<10 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      HistMod_t::MakeLine(HistMod_t::GetColor(mctrue-9), 1.0)
                                                      );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 10)
                this->GetHist(4).Fill(f);
        }
    }
};

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct OmegaHist_t {

    using Tree_t = physics::OmegaEtaG2::OmegaTree_t;

    struct Fill_t {
        const Tree_t& Tree;
        Fill_t(const Tree_t& t) : Tree(t) {}
        double TaggW() const {
            return Tree.TaggW;
        }
    };

    TH1D* h_KinFitChi2;
    TH1D* h_gggIM;
    TH1D* h_ggIM;
    TH1D* h_mm;
    TH1D* h_bachelorE;
    TH2D* h_p_Theta_E;
    TH2D* h_mm_gggIM;

    const BinSettings Ebins = BinSettings(1600,   0, 1600);

    const BinSettings Chi2bins = BinSettings (250, 0,   25);
    const BinSettings probbins = BinSettings (250, 0,   1);

    const BinSettings IMbins = BinSettings(1600,   0, 1600);
    const BinSettings MMbins = BinSettings(1600, 400, 2000);

    const BinSettings MMgggIMbins_X = BinSettings(600, 0, 1200);
    const BinSettings MMgggIMbins_Y = BinSettings(750, 500, 2000);

    const BinSettings pThetaBins = BinSettings(500, 0, 50);
    const BinSettings pEbins = BinSettings(1000,   0, 1000);

    OmegaHist_t(HistogramFactory HistFac) {
        h_KinFitChi2 = HistFac.makeTH1D("KinFitChi2", "#chi^{2}",             "", BinSettings(200,0,100), "h_KinFitChi2");
        h_gggIM      = HistFac.makeTH1D("3#gamma IM", "3#gamma IM [MeV]",     "", IMbins,                 "h_ggg_IM");
        h_ggIM       = HistFac.makeTH1D("2#gamma sub-IM", "2#gamma IM [MeV]", "", IMbins,                 "h_gg_IM");
        h_mm         = HistFac.makeTH1D("Missing Mass",   "MM [MeV]",         "", IMbins,                 "h_mm");
        h_bachelorE  = HistFac.makeTH1D("Bachelor E",     "E [MeV]",          "", IMbins,                 "h_bachelorE");
        h_p_Theta_E  = HistFac.makeTH2D("Proton #theta vs. E_{k}", "E_{k} [MeV]", "#theta [#circ]", pEbins, pThetaBins, "h_p_theta_E");
        h_mm_gggIM   = HistFac.makeTH2D("Missing Mass / 3#gamma IM", "3#gamma IM [MeV]", "MM [MeV]", IMbins, MMbins, "h_mm_gggIM");
    }

    void Fill(const Fill_t& f) const {

        h_KinFitChi2->Fill(f.Tree.KinFitChi2, f.TaggW());

        h_gggIM->Fill(f.Tree.ggg().M(), f.TaggW());

        for(const auto& v : f.Tree.ggIM())
            h_ggIM->Fill(v, f.TaggW());

        h_mm->Fill(f.Tree.mm().M(), f.TaggW());

        for(const auto& v : f.Tree.BachelorE())
            h_bachelorE->Fill(v, f.TaggW());

        h_p_Theta_E->Fill(f.Tree.ggg().M(), radian_to_degree(f.Tree.p().Theta()));

        h_mm_gggIM->Fill(f.Tree.ggg().M(), f.Tree.mm().M(), f.TaggW());

    }

    std::vector<TH1*> GetHists() const {
        return {h_KinFitChi2,h_gggIM,h_ggIM,h_mm,h_bachelorE,h_p_Theta_E,h_mm_gggIM};
    }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                                 {"KinFitChi2<10", [] (const Fill_t& f) { return f.Tree.KinFitChi2< 10; } },
                                 {"mm cut",        [] (const Fill_t& f) { return f.Tree.mm().M()<1100 && f.Tree.mm().M() > 800; } },
                             });
        return cuts;
    }

};



int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("plot", ' ', "0.1");
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","input");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>("m","maxevents","Process only max events",false,"maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("","tree","Tree name",false,"Fitted/SigAll","treename");

    cmd.parse(argc, argv);

    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-b");
    }
    TRint app("EtapOmegaG_plot",&fake_argc,fake_argv,nullptr,0,true);
    auto oldsig = app.GetSignalHandler();
    oldsig->Remove();
    auto mysig = new MyTInterruptHandler();
    mysig->Add();
    gSystem->AddSignalHandler(mysig);

    WrapTFileInput input(cmd_input->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        if(wraptree->Matches(t)) {
            wraptree->LinkBranches(t);
            return true;
        }
        return false;
    };


    OmegaHist_t::Tree_t tree;

    if(!link_branches("OmegaEtaG2/tree", addressof(tree), -1)) {
        LOG(ERROR) << "Cannot link branches of tree";
        return 1;
    }

    const auto entries = tree.Tree->GetEntries();

     unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                    WrapTFileOutput::mode_t::recreate,
                                                     true); // cd into masterFile upon creation
    }


    HistogramFactory HistFac("OmegaEtaG2");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");

    auto cuttree = cuttree::Make<MCTrue_Splitter<OmegaHist_t>>(HistFac,
                                              sanitized_treename,
                                              OmegaHist_t::GetCuts()
                                              );

    LOG(INFO) << "Tree entries=" << entries;

    auto max_entries = entries;
    if(cmd_maxevents->isSet() && cmd_maxevents->getValue().back()<entries) {
        max_entries = cmd_maxevents->getValue().back();
        LOG(INFO) << "Running until " << max_entries;
    }

    for(long long entry=0;entry<max_entries;entry++) {
        if(interrupt)
            break;

        tree.Tree->GetEntry(entry);
        cuttree::Fill<MCTrue_Splitter<OmegaHist_t>>(cuttree, {tree});

        if(entry % 100000 == 0)
            LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

            mysig->Remove();
            oldsig->Add();
            gSystem->AddSignalHandler(oldsig);
            delete mysig;

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
