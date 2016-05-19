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
#include "base/vec/vec2.h"

#include <list>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TRint.h"

#include "TStyle.h"

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

        using histstyle::Mod_t;

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
                                                       Mod_t::MakeLine(histstyle::color_t::Get(mctrue-10), 1, kGray+1)
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

        int iBestIndex() const {
            if(Tree.bestHyp == 2) {
                return Tree.iBestEta;
            } else if(Tree.bestHyp == 1) {
                return Tree.iBestPi0;
            } else
                return -1;
        }

        double BestBachelorE() const {
            return iBestIndex() != -1 ? Tree.BachelorE_fitted().at(iBestIndex()) : NaN;
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

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings Ebins    = Bins(1600, 0, 1600);

    const BinSettings Chi2Bins = BinSettings(250, 0,   25);
    const BinSettings probbins = BinSettings(250, 0,   1);

    const BinSettings IMbins        = Bins(1600,   0, 1600);
    const BinSettings gggIMbins     = Bins( 300,   650, 950);
    const BinSettings MMbins        = Bins(1600, 400, 2000);

    const BinSettings MMgggIMbins_X = Bins(600, 0, 1200);
    const BinSettings MMgggIMbins_Y = Bins(750, 500, 2000);

    const BinSettings pThetaBins = Bins( 125, 0,   50);
    const BinSettings pEbins     = Bins(250,  0, 1000);
    const BinSettings PSAABins   = Bins(  60, 20,   60);
    const BinSettings PSARBins   = Bins( 100,  0,  450);
    const BinSettings TaggChBins = BinSettings(47);

    const BinSettings TaggTime   = BinSettings(200, -25, 25);
    const BinSettings CoplBins   = Bins(300, 0, 30.0);

    const BinSettings dalitzBins = Bins(200, -0.4, 0.4);

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
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.TaggW());
        });

        AddTH1("2#gamma sub-IM",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM",
               [] (TH1D* h, const Fill_t& f) {

            for(const auto& v : f.Tree.ggIM_fitted())
                h->Fill(v, f.TaggW());
        });

        AddTH1("Bachelor Photon Energy",  "E [MeV]",     "",       IMbins,     "bachelorE",
               [] (TH1D* h, const Fill_t& f) {

            for(const auto& v : f.Tree.BachelorE_fitted())
                h->Fill(v, f.TaggW());
        });

        AddTH1("Bachelor Photon Energy|Best Hyp",  "E [MeV]",     "",       IMbins,     "bachelorE_bestHyp",
               [] (TH1D* h, const Fill_t& f) {

            if( f.iBestIndex() != -1 ) {
                h->Fill(f.Tree.BachelorE_fitted().at(f.iBestIndex()), f.TaggW());
            }
        });

        AddTH1("2#gamma sub-IM|Best Hyp",  "2#gamma IM [MeV]",     "",       IMbins,     "gg_IM_bestHyp",
               [] (TH1D* h, const Fill_t& f) {

            if( f.iBestIndex() != -1 ) {
                h->Fill(f.Tree.ggIM().at(f.iBestIndex()), f.TaggW());
            }
        });


        AddTH1("Missing Mass",      "MM [MeV]",     "",       MMbins,     "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });


//        AddTH1("Tagger Channels", "Channel",              "# hits", TaggChBins, "TaggCh",
//               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggCh, f.TaggW());
//        });

//        AddTH1("Coplanarity Angle", "Coplanarity angle [#circ]", "", CoplBins, "CoplAngle",
//               [] (TH1D* h, const Fill_t& f) { h->Fill(radian_to_degree(f.Tree.copl_angle()), f.TaggW());
//        });

        AddTH1("Tagger Time - CB Average Time", "t [ns]", "",       TaggTime,   "TaggTime",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
        });

//        AddTH2("Proton #theta vs. E_{k}", "E_{k} [MeV]", "#theta [#circ]",  pEbins,   pThetaBins, "p_theta_E",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.p_fitted().E() - ParticleTypeDatabase::Proton.Mass(), radian_to_degree(f.Tree.p_fitted().Theta()), f.TaggW());
//        });

//        AddTH2("Missing Mass / 3#gamma IM", "3#gamma IM [MeV]", "MM [MeV]", IMbins,   MMbins,     "mm_gggIM",
//               [] (TH2D* h, const Fill_t& f) { h->Fill(f.Tree.ggg_fitted().M(), f.Tree.mm().M(), f.TaggW());
//        });

//        AddTH2("Proton PSA", "PSA Angle [#circ]", "PSA Radius",             PSAABins, PSARBins,   "p_PSA",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.p_PSA_Angle, f.Tree.p_PSA_Radius, f.TaggW());
//        });

//        AddTH2("3#gamma IM vs 2#gamma IM", "3#gamma IM [MeV]", "max(2#gamma IM) [MeV]", IMbins, IMbins, "ggg_max_gg",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.ggg_fitted().M(), max(f.Tree.ggIM_fitted()), f.TaggW());
//        });

        AddTH1("#pi^{0} Hyp: prob", "prob_{#pi^{0}}","",probbins, "pi0hyp_prob",
               [] (TH1D* h, const Fill_t& f) {
            const auto& i = f.Tree.iBestPi0;
            if(i >= 0)
                h->Fill(f.Tree.pi0prob().at(size_t(i)), f.TaggW());
        });

        AddTH1("#eta Hyp: prob", "#chi^{2}_{#eta}","",probbins, "etahyp_prob",
               [] (TH1D* h, const Fill_t& f) {
            const auto& i = f.Tree.iBestEta;
            if(i >= 0)
                h->Fill(f.Tree.etaprob().at(size_t(i)), f.TaggW());
        });

        AddTH1("#eta Hyp: #omega IM", "m(#omega_{#eta})","",IMbins, "etahyp_omega",
               [] (TH1D* h, const Fill_t& f) {
            const auto& i = f.Tree.iBestEta;
            if(i >= 0)
                h->Fill(f.Tree.eta_omega_im().at(size_t(i)), f.TaggW());
        });

        AddTH1("#pi^{0} Hyp: #omega IM", "m(#omega_{#pi^{0}}})","",IMbins, "pi0hyp_omega",
               [] (TH1D* h, const Fill_t& f) {
            const auto& i = f.Tree.iBestEta;
            if(i >= 0)
                h->Fill(f.Tree.pi0_omega_im().at(size_t(i)), f.TaggW());
        });

        AddTH1("Extra Gamma makes #pi^{0}", "m(2#gamma)", "", IMbins, "extaGammaPi0",
                [] (TH1D* h, const Fill_t& f) {

            h->Fill(f.Tree.bachelor_extra().at(f.iBestIndex()).M(), f.TaggW());
        });

//        AddTH1("verify", "E", "", IMbins, "verify",
//                [] (TH1D* h, const Fill_t& f) {
//            const vector<int> index = {0,1,2};
//            auto boost = - f.Tree.ggg().BoostVector();
//            auto ba = f.Tree.photons().at(index.at(f.iBestIndex()));
//            ba.Boost(boost);

//            h->Fill(ba.E(), f.TaggW());
//        });

        AddTH2("Dalitz","X","Y", dalitzBins, dalitzBins, "dalitz",
               [] (TH2D* h, const Fill_t& f) {

            OmegaDalitzPlot p(f.Tree.photons_fitted(), f.Tree.ggg_fitted());
            do {
                h->Fill(p.var.x, p.var.y);
            } while (p.Next());

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
                                 {"Prob>0.002+mm", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.002 && f.Tree.mm().M()<1100 && f.Tree.mm().M() > 780; } }
                             });

        auto etaHypCut = [] (const Fill_t& f) {

            if(f.Tree.bestHyp != 2)
                return false;

            const auto& etaprob  = f.Tree.etaprob()[f.Tree.iBestEta];
            const auto& iBestPi0 = f.Tree.iBestPi0;
            const auto& pi0prob  = f.Tree.pi0prob()[iBestPi0];

            return etaprob > 0.03 && (iBestPi0==-1 || pi0prob < 0.03);
        };

        auto etaHypCut2 = [] (const Fill_t& f) {

            if(f.Tree.bestHyp != 2)
                return false;

            const auto& etaprob  = f.Tree.etaprob()[f.Tree.iBestEta];
            const auto& iBestPi0 = f.Tree.iBestPi0;
            const auto& pi0prob  = f.Tree.pi0prob()[iBestPi0];

            return etaprob > 0.4 && (iBestPi0==-1 || pi0prob < 0.03);
        };

        auto pi0HypCut = [] (const Fill_t& f) {
            if(f.Tree.bestHyp != 1)
                return false;

            const auto& iBestPi0 = f.Tree.iBestPi0;
            const auto& pi0prob  = f.Tree.pi0prob()[iBestPi0];

            return pi0prob > 0.03;
        };

        auto etaBachelorCut = [] (const Fill_t& f) {
            return interval<double>::CenterWidth(200,40).Contains(f.BestBachelorE());
        };

        auto pi0BachelorCut = [] (const Fill_t& f) {
            return interval<double>::CenterWidth(380,40).Contains(f.BestBachelorE());
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"m(3#gamma) cut",        [] (const Fill_t& f) { return f.Tree.ggg_fitted().M()<900 && f.Tree.ggg_fitted().M() > 700; } },
                              {"etaHyp",               etaHypCut},
                              {"etaHyp2",              etaHypCut2},
                              {"pi0Hyp",               pi0HypCut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"eta-bachelorCut",      etaBachelorCut},
                              {"pi0-bachelorCut",      pi0BachelorCut}
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

    auto cmd_tree = cmd.add<TCLAP::ValueArg<string>>("","tree","Tree name",false,"T1","treename");
    auto cmd_pres = cmd.add<TCLAP::SwitchArg>("p","", "Presentation Mode", false);
    auto cmd_binscale = cmd.add<TCLAP::ValueArg<double>>("B","bin-scale","Bin Scaleing", false, 1.0, "bins");

    cmd.parse(argc, argv);

    if(cmd_pres->isSet()) {
        gStyle->SetLabelSize(.05f, "XYZ");
        gStyle->SetTextSize(.05f);
        gStyle->SetCanvasBorderSize(0);
    }

    if(cmd_binscale->isSet()) {
        binScale = cmd_binscale->getValue();
    }


    WrapTFileInput input(cmd_input->getValue());

    auto link_branches = [&input] (const string treename, WrapTTree* wraptree, long long expected_entries) {
        TTree* t;
        if(!input.GetObject(treename,t))
            throw runtime_error("Cannot find tree "+treename+" in input file");
        if(expected_entries>=0 && t->GetEntries() != expected_entries)
            throw runtime_error("Tree "+treename+" does not have entries=="+to_string(expected_entries));
        if(wraptree->Matches(t,false)) {
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

    auto signal_hists = cuttree::Make<MCTrue_Splitter<OmegaHist_t>>(HistFac,
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
        cuttree::Fill<MCTrue_Splitter<OmegaHist_t>>(signal_hists, {tree});

        if(entry % 100000 == 0)
            LOG(INFO) << "Processed " << 100.0*entry/entries << " %";
    }

    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {

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
