/**
  * @file OmegaEtaG2_plot.cc
  * @brief Plotting tool for TTrees from the ant::analysis::physics::OmegaEtaG2 physics class
  */

#include "base/Logger.h"

#include "analysis/plot/CutTree.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/production/triplePi0.h"

#include "analysis/utils/particle_tools.h"

#include "base/CmdLine.h"
#include "base/interval.h"
#include "base/printable.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/math.h"
#include "base/vec/vec2.h"

#include <list>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TRint.h"

#include "TStyle.h"
#include "TCutG.h"

#include "TChain.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

volatile static bool interrupt = false;
static double binScale=1.0;

class OmegaDalitzPlot {

public:

    static double x(double T1, double T2, double T3) noexcept {
        return (T2 - T1) / (sqrt(2) * (T1+T2+T3));
    }

    static double y(double T1, double T2, double T3) noexcept {
        return T3 / (T1+T2+T3) - 1.0/3.0;
    }

    static vec2 xy(double T1, double T2, double T3) noexcept {
        return vec2(x(T1,T2,T3),y(T1,T2,T3));
    }

    static std::vector<double> getDalitzT(const std::vector<TLorentzVector>& photons, const TLorentzVector& omega) {

        const auto boost = -omega.BoostVector();

        vector<double> T;
        T.reserve(photons.size());
        for( const auto& g : photons) {
            TLorentzVector lv = g;
            lv.Boost(boost);
            T.push_back(lv.E());
        }

        sort(T.begin(),T.end());

        return T;
    }

protected:
    vec2 calc() {
        return xy(T.at(0), T.at(1), T.at(2));
    }

    std::vector<double> T;

public:
    vec2 var;

    OmegaDalitzPlot(const std::vector<TLorentzVector>& photons, const TLorentzVector& omega):
        T(getDalitzT(photons,omega)),
        var(calc())
    {}

    bool Next() noexcept {
        const auto r = std::next_permutation(T.begin(), T.end());
        if(r)
            var = calc();
        return r;
    }

};

static string data_name = "Data";

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;



    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)

    {
        using histstyle::Mod_t;

        // TODO: derive this from channel map
        this->GetHist(0, data_name, Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2));
        this->GetHist(2, "MainBkg",  Mod_t::MakeLine(kGreen, 2));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = f.Tree.MCTrue;

        using histstyle::Mod_t;

        auto get_bkg_name = [] (const unsigned mctrue) {
            return physics::triplePi0::getOtherChannelNames(mctrue); //(int(mctrue));
        };

        using histstyle::Mod_t;

        const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                       get_bkg_name(mctrue),
                                                       Mod_t::MakeLine(histstyle::color_t::Get(mctrue-10), 1, kGray+1)
                                                       );


        hist.Fill(f);

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 9 || mctrue == 2)
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

double maxIM(const std::vector<TLorentzVector>& data) {
    return max_element(data.cbegin(), data.cend(), [] (const TLorentzVector& v1, const TLorentzVector& v2) { return v1.M() < v2.M(); })->M();
}

double maxE(const std::vector<TLorentzVector>& data) {
    return max_element(data.cbegin(), data.cend(), [] (const TLorentzVector& v1, const TLorentzVector& v2) { return v1.E() < v2.E(); })->M();
}

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct TriplePi0Hist_t {

    using Tree_t = physics::triplePi0::PionProdTree;

    struct Fill_t {
        const Tree_t& Tree;

        Fill_t(const Tree_t& t) : Tree(t) {}

        double TaggW() const {
            return Tree.Tagg_W;
        }

//        int iBestIndex() const {
//            if(Tree.bestHyp == 2) {
//                return Tree.iBestEta;
//            } else if(Tree.bestHyp == 1) {
//                return Tree.iBestPi0;
//            } else
//                return -1;
//        }

//        double BestBachelorE() const {
//            return iBestIndex() != -1 ? Tree.BachelorE_fitted().at(size_t(iBestIndex())) : NaN;
//        }

//        int BachelorIndex() const {
//            return iBestIndex();
//        }
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

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings Ebins    = Bins(1000, 0, 1000);

    const BinSettings Chi2Bins = BinSettings(250, 0,   25);
    const BinSettings probbins = BinSettings(250, 0,   1);

    const BinSettings IMbins        = Bins(1000,  200, 1100);
    const BinSettings IMProtonBins  = Bins(1000,  600, 1200);
    const BinSettings IM2g          = Bins(1000,    0,  360);

    const BinSettings pThetaBins = Bins( 200,  0,   80);
    const BinSettings pEbins     = Bins( 350,  0, 1200);

    HistogramFactory HistFac;

    void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, fillfunc_t<TH1D> f) {
        h1.emplace_back(HistFiller_t<TH1D>(
                            HistFac.makeTH1D(title, xlabel, ylabel, bins, name),f));
    }

    void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, fillfunc_t<TH2D> f) {
        h2.emplace_back(HistFiller_t<TH2D>(
                            HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name),f));
    }

    TriplePi0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf) {



        AddTH1("KinFit Probability",      "probability",             "",       probbins,   "KinFitProb",
               [] (TH1D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.EMB_prob, f.TaggW());
        });

        AddTH1("6#gamma IM","6#gamma IM [MeV]", "", IMbins,"IM_6g",
               [] (TH1D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.IM6g, f.TaggW());
        });
        AddTH1("6#gamma IM fitted","6#gamma IM [MeV]", "", IMbins,"IM_6g_fit",
               [] (TH1D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.EMB_IM6g, f.TaggW());
        });

        AddTH1("MM proton","MM_{proton} [MeV]", "", IMProtonBins, "IM_p",
               [] (TH1D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.proton_MM().M(), f.TaggW());
        });

        AddTH1("MM proton fitted","MM_{proton} [MeV]", "", IMProtonBins, "IM_p_fit",
               [] (TH1D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.EMB_proton_MM().M(),f.TaggW());
        });

        AddTH1("MM pions", "IM_{2#gamma} [MeV]","", IM2g,"IM_pions",
               [] (TH1D* h, const Fill_t& f)
        {
            for ( const auto& pion: f.Tree.SIG_pions())
                h->Fill(pion.M(),f.TaggW());
        });

        AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE",
               [] (TH2D* h, const Fill_t& f)
        {
            h->Fill(f.Tree.EMB_proton().E() - 938.3, std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
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

    static TCutG* makeDalitzCut() {
        TCutG* c = new TCutG("DalitzCut", 3);
        c->SetPoint(0, 0.0,  0.2);
        c->SetPoint(1, -.22, -.11);
        c->SetPoint(2,  .22, -.11);
        return c;
    }

    static TCutG* dalitzCut;

    struct TreeCuts {

        static bool KinFitProb_MM(const Fill_t& f) noexcept {
            return     f.Tree.EMB_prob >  0.05;
        }
    };

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;


        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"EMB_prob > 0.1",  [](const Fill_t& f)
                                           {
                                               return f.Tree.EMB_prob > 0.1;
                                           }
                              }
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"SIG_prob > 0.1",  [](const Fill_t& f)
                                           {
                                               return f.Tree.SIG_prob > 0.1;
                                           }
                              }
                          });
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"BKG_prob < 0.1",  [](const Fill_t& f)
                                           {
                                               return f.Tree.BKG_prob < 0.1;
                                           }
                              }
                          });


        return cuts;
    }

};



int main(int argc, char** argv) {
    SetupLogger();

    signal(SIGINT, [] (int) { interrupt = true; } );

    TCLAP::CmdLine cmd("plot", ' ', "0.1");

    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch",    "Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_maxevents = cmd.add<TCLAP::MultiArg<int>>( "m", "maxevents","Process only max events",                     false, "maxevents");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>( "o", "output",    "Output file",                                false, "",            "filename");
    auto cmd_inputFile = cmd.add<TCLAP::ValueArg<string>>( "i", "input",    "Input file",                                false, "",            "filename");

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("","tree","Tree name",false,"T1","treename");
    auto cmd_pres = cmd.add<TCLAP::SwitchArg>("p","", "Presentation Mode", false);
    auto cmd_binscale = cmd.add<TCLAP::ValueArg<double>>("B","bin-scale","Bin Scaleing", false, 1.0, "bins");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_dataName = cmd.add<TCLAP::ValueArg<string>>("","DataName","Name of the data set",false,"Data","dataname");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    if(cmd_pres->isSet()) {
        gStyle->SetLabelSize(.05f, "XYZ");
        gStyle->SetTextSize(.05f);
        gStyle->SetCanvasBorderSize(0);
    }

    if(cmd_binscale->isSet()) {
        binScale = cmd_binscale->getValue();
    }

    if(cmd_dataName->isSet()) {
        data_name = cmd_dataName->getValue();
    }


    const string treename = physics::triplePi0::treeAccessName();

    TriplePi0Hist_t::Tree_t tree;

    WrapTFileInput infile(cmd_inputFile->getValue());

    TTree* t;

    infile.GetObject(treename, t);


    if(!tree.Matches(t,false))
        throw runtime_error("Tree branches don't match");

    tree.LinkBranches(t);

    const auto entries = tree.Tree->GetEntries();

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(),
                                                           WrapTFileOutput::mode_t::recreate,
                                                           true); // cd into masterFile upon creation
    }


    HistogramFactory HistFac("OmegaEtaG2");

    const auto& sanitized_treename = std_ext::replace_str(cmd_tree->getValue(),"/","_");

    auto signal_hists = cuttree::Make<MCTrue_Splitter<TriplePi0Hist_t>>(HistFac,
                                                                    sanitized_treename,
                                                                    TriplePi0Hist_t::GetCuts()
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
        cuttree::Fill<MCTrue_Splitter<TriplePi0Hist_t>>(signal_hists, {tree});

        if(entry % 100000 == 0)
            LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
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

TCutG* TriplePi0Hist_t::dalitzCut = TriplePi0Hist_t::makeDalitzCut();
