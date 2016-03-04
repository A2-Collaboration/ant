/**
  * @file OmegaEtaG2_plot.cc
  * @brief Plotting tool for TTrees from the ant::analysis::physics::OmegaEtaG2 physics class
  */

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

#include <list>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TRint.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile static bool interrupt = false;

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
        using histstyle::Mod_t;

        // TODO: derive this from channel map
        this->GetHist(0, "Data", Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2));
        this->GetHist(2, "Ref",  Mod_t::MakeLine(kGreen, 2));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = unsigned(f.Tree.Channel);

        auto get_bkg_name = [] (const unsigned mctrue) {
            const auto entry = physics::OmegaEtaG2::reaction_channels.channels.find(int(mctrue));

            if(entry!=physics::OmegaEtaG2::reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        using histstyle::Mod_t;
        const Hist_t& hist = mctrue<10 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      Mod_t::MakeLine(histstyle::color_t::Get(mctrue-9), 1)
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

bool Contains(const interval<double>& i, const std::vector<double>& d) {
    for(const auto& v : d)
        if(i.Contains(v))
            return true;

    return false;
}


double max(const std::vector<double>& data) {
    return *max_element(data.cbegin(), data.cend());
}


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

    template <typename Hist>
    using fillfunc_t = std::function<void(Hist*, const Fill_t&)>;

    template <typename Hist>
    struct HistFiller_t {
        fillfunc_t<Hist> func;
        Hist* h;
        HistFiller_t(Hist* hist, fillfunc_t<Hist> f): func(f), h(hist) {}
        void Fill(const Fill_t& data) const {
            func(h, data);
        }
    };

    template <typename Hist>
    struct HistMgr : std::list<HistFiller_t<Hist>> {

        using list<HistFiller_t<Hist>>::list;

        void Fill(const Fill_t& data) const {
            for(auto& h : *this) {
                h.Fill(data);
            }
        }
    };

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings Ebins = BinSettings(1600,   0, 1600);

    const BinSettings Chi2bins = BinSettings (2500, 0,   25);
    const BinSettings probbins = BinSettings (250, 0,   1);

    const BinSettings IMbins = BinSettings(1600,   0, 1600);
    const BinSettings MMbins = BinSettings(1600, 400, 2000);

    const BinSettings MMgggIMbins_X = BinSettings(600, 0, 1200);
    const BinSettings MMgggIMbins_Y = BinSettings(750, 500, 2000);

    const BinSettings pThetaBins = BinSettings( 500,  0,   50);
    const BinSettings pEbins     = BinSettings(1000,  0, 1000);
    const BinSettings PSAABins   = BinSettings(  60, 20,   60);
    const BinSettings PSARBins   = BinSettings( 100,  0,  450);
    const BinSettings TaggChBins = BinSettings(47);
    const BinSettings Chi2Bins   = BinSettings(100,0,50);
    const BinSettings TaggTime   = BinSettings(200, -25, 25);
    const BinSettings CoplBins   = BinSettings(300, 0, 30.0);

    HistogramFactory HistFac;

    void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, fillfunc_t<TH1D> f) {
        h1.emplace_back(HistFiller_t<TH1D>(
                                    HistFac.makeTH1D(title, xlabel, ylabel, bins, name),f));
    }

    void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, fillfunc_t<TH2D> f) {
        h2.emplace_back(HistFiller_t<TH2D>(
                                    HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name),f));
    }

    OmegaHist_t(const HistogramFactory& hf): HistFac(hf) {

        AddTH1("KinFitChi2",      "#chi^{2}",             "",       Chi2Bins,   "KinFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.KinFitChi2, f.TaggW());
        });

        AddTH1("3#gamma IM",      "3#gamma IM [MeV]",     "",       IMbins,     "ggg_IM",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg().M(), f.TaggW());
        });

        AddTH1("2#gamma sub-IM",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM",
               [] (TH1D* h, const Fill_t& f) {

            for(const auto& v : f.Tree.ggIM())
                h->Fill(v, f.TaggW());
        });

        AddTH1("Bachelor Photon Energy",  "E [MeV]",     "",       IMbins,     "bachelorE",
               [] (TH1D* h, const Fill_t& f) {

            for(const auto& v : f.Tree.BachelorE())
                h->Fill(v, f.TaggW());
        });

        AddTH1("Missing Mass",      "MM [MeV]",     "",       MMbins,     "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });


        AddTH1("Tagger Channels", "Channel",              "# hits", TaggChBins, "TaggCh",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggCh, f.TaggW());
        });

        AddTH1("Coplanarity Angle", "Coplanarity angle [#circ]", "", CoplBins, "CoplAngle",
               [] (TH1D* h, const Fill_t& f) { h->Fill(radian_to_degree(f.Tree.copl_angle()), f.TaggW());
        });

        AddTH1("Tagger Time - CB Average Time", "t [ns]", "",       TaggTime,   "TaggTime",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
        });

        AddTH2("Proton #theta vs. E_{k}", "E_{k} [MeV]", "#theta [#circ]",  pEbins,   pThetaBins, "p_theta_E",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), radian_to_degree(f.Tree.p_fitted().Theta()));
        });

        AddTH2("Missing Mass / 3#gamma IM", "3#gamma IM [MeV]", "MM [MeV]", IMbins,   MMbins,     "mm_gggIM",
               [] (TH2D* h, const Fill_t& f) { h->Fill(f.Tree.ggg().M(), f.Tree.mm().M(), f.TaggW());
        });

        AddTH2("Proton PSA", "PSA Angle [#circ]", "PSA Radius",             PSAABins, PSARBins,   "p_PSA",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.p_PSA_Angle, f.Tree.p_PSA_Radius, f.TaggW());
        });

        AddTH2("3#gamma IM vs 2#gamma IM", "3#gamma IM [MeV]", "max(2#gamma IM) [MeV]", IMbins, IMbins, "ggg_max_gg",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.ggg().M(), max(f.Tree.ggIM()), f.TaggW());
        });

    }

    void Fill(const Fill_t& f) const {

        h1.Fill(f);
        h2.Fill(f);

    }

    std::vector<TH1*> GetHists() const {
        vector<TH1*> v;
        v.reserve(h1.size()+h2.size());
        for(auto& e : h1) {
            v.emplace_back(e.h);
        }
        for(auto& e: h2) {
            v.emplace_back(e.h);
        }
        return v;
    }

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                                 {"KinFitChi2<5 ", [] (const Fill_t& f) { return f.Tree.KinFitChi2<5; } }
                             });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"mm cut",        [] (const Fill_t& f) { return f.Tree.mm().M()<1100 && f.Tree.mm().M() > 780; } }
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"gggIM cut",        [] (const Fill_t& f) { return f.Tree.ggg().M()<900 && f.Tree.ggg().M() > 700; } },
                              {"eta",           [] (const Fill_t& f) { return Contains( {530.0, 580.0}, f.Tree.ggIM()); } },
                              {"pi0",           [] (const Fill_t& f) { return Contains( {125.0, 145.0}, f.Tree.ggIM()); } }

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

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("","tree","Tree name",false,"T1","treename");

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

    canvas c("Steps");
    hstack s;
    for(const auto& c : physics::OmegaEtaG2::reaction_channels.channels) {
        TH1D* h = nullptr;
        cout << "OmegaEtaG2/setps_"+to_string(c.first) << endl;
        input.GetObject("OmegaEtaG2/steps_"+to_string(c.first), h);
        if(h){
            s << h;
            cout << h->GetName() << endl;
        }
    }
    c << &s << endc;

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
