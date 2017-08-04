#include "MesonDalitzDecays.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "base/Logger.h"
#include "base/interval.h"
#include "base/std_ext/vector.h"

#include "expconfig/ExpConfig.h"


using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;


template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    const decltype(physics::MesonDalitzDecays::makeChannels()) channels;

    constexpr static Color_t bkg_color = kGray+1;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
        channels(physics::MesonDalitzDecays::makeChannels())
    {
        using histstyle::Mod_t;

        this->GetHist(0, "Data",   Mod_t::MakeDataPoints(kBlack));
        this->GetHist(1, "Signal", Mod_t::MakeLine(kRed, 2));
        // mctrue is never >= 3 (and < 9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(2, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(3, "Bkg_MC", Mod_t::MakeFill(bkg_color, -1));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = unsigned(f.Tree.channel);

        using histstyle::Mod_t;

        auto get_bkg_name = [] (const unsigned mctrue) {
            const auto entry = physics::MesonDalitzDecays::reaction_channels.channels.find(int(mctrue));

            if (entry != physics::MesonDalitzDecays::reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        auto get_color = [] (const unsigned mctrue) -> short {
            const auto entry = physics::MesonDalitzDecays::reaction_channels.channels.find(int(mctrue));

            if (entry != physics::MesonDalitzDecays::reaction_channels.channels.end())
                return entry->second.color;

            return histstyle::color_t::GetLight(mctrue-10);
        };

        const Hist_t& hist = mctrue < 10 ? this->GetHist(mctrue)
                                         : this->GetHist(mctrue,
                                                         get_bkg_name(mctrue),
                                                         Mod_t::MakeLine(get_color(mctrue), 1, bkg_color)
                                                         );

        hist.Fill(f);

        // handle MC_all and MC_bkg
        if (mctrue > 0) {
            this->GetHist(2).Fill(f);
            if (mctrue >= 10)
                this->GetHist(3).Fill(f);
        }
    }
};

bool Contains(const interval<double>& i, const std::vector<double>& d)
{
    for (const auto& v : d)
        if (i.Contains(v))
            return true;

    return false;
}

double max(const std::vector<double>& data)
{
    return *max_element(data.cbegin(), data.cend());
}

double im_ee(vector<double> vetoE, vector<TLorentzVector> photons)
{
    const auto leptons = std_ext::get_sorted_indices(vetoE);

    return (photons.at(leptons[0]) + photons.at(leptons[1])).M();
}


// define the structs containing the histograms and the cuts
struct Hist_t {

    using Tree_t = physics::MesonDalitzDecays::Tree_t;

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
            for (auto& h : *this)
                h.Fill(data);
        }
    };

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    const BinSettings pid_channels = BinSettings(
                ExpConfig::Setup::GetDetector(Detector_t::Type_t::PID)->GetNChannels());

    const BinSettings Ebins      = BinSettings(1200, 0, 1200);

    const BinSettings Chi2Bins   = BinSettings(250, 0, 25);
    const BinSettings probbins   = BinSettings(250, 0,  1);

    const BinSettings IMbins     = BinSettings(1200,   0, 1200);
    const BinSettings MMbins     = BinSettings(1200, 400, 1600);

    const BinSettings pThetaBins = BinSettings( 125,  0,   50);
    const BinSettings pEbins     = BinSettings( 250,  0, 1000);
    const BinSettings PSAABins   = BinSettings(  60, 20,   60);
    const BinSettings PSARBins   = BinSettings( 100,  0,  450);
    const BinSettings TaggChBins = BinSettings(47);

    const BinSettings TaggTime   = BinSettings(240, -30, 30);
    const BinSettings CoplBins   = BinSettings(300, 0, 30);

    const BinSettings zVertex    = BinSettings(100, -15, 15);

    const BinSettings vetoEbins  = BinSettings(200, 0, 10);

    HistogramFactory HistFac;

    void AddTH1(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &bins, const string &name, fillfunc_t<TH1D> f)
    {
        h1.emplace_back(HistFiller_t<TH1D>(
                            HistFac.makeTH1D(title, xlabel, ylabel, bins, name), f));
    }

    void AddTH2(const string &title, const string &xlabel, const string &ylabel,
                const BinSettings &xbins, const BinSettings& ybins,
                const string &name, fillfunc_t<TH2D> f)
    {
        h2.emplace_back(HistFiller_t<TH2D>(
                            HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name), f));
    }

    Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf) {

        AddTH1("KinFitChi2", "#chi^{2}", "#", Chi2Bins, "KinFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_chi2, f.TaggW());
        });

        AddTH1("TreeFitChi2", "#chi^{2}", "#", Chi2Bins, "TreeFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_chi2, f.TaggW());
        });

        AddTH1("KinFitProb", "probability", "#", probbins, "KinFitProb",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_probability, f.TaggW());
        });

        AddTH1("TreeFitProb", "probability", "#", probbins, "TreeFitProb",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_probability, f.TaggW());
        });

        AddTH1("3 photon IM", "3#gamma IM [MeV]", "#", IMbins, "etaIM",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.eta().M(), f.TaggW());
        });

        AddTH1("3 photon IM kinfitted",  "3#gamma IM fit [MeV]", "#", IMbins, "etaIM_kinfitted",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.eta_kinfit().M(), f.TaggW());
        });

//        AddTH1("3 photon IM treefitted",  "3#gamma IM fit [MeV]", "#", IMbins, "etaIM_treefitted",
//               [] (TH1D* h, const Fill_t& f) {
//            h->Fill(f.Tree.eta_treefit().M(), f.TaggW());
//        });

        AddTH1("Missing Mass", "MM [MeV]", "", MMbins, "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });

        AddTH1("Z Vertex Kinfit", "z [cm]", "#", zVertex, "v_z_kinfit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex Treefit", "z [cm]", "#", zVertex, "v_z_treefit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_ZVertex, f.TaggW());
        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.eta().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()), f.TaggW());
//        });

        AddTH2("IM(e+e-) vs. IM(e+e-g) fit", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.eta_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()), f.TaggW());
        });

        AddTH2("PID 2 charged 1 neutral", "PID Energy [MeV]", "PID Channel", vetoEbins, pid_channels, "eegPID",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_vetoE().at(i), f.Tree.photons_vetoChannel().at(i), f.TaggW());
        });

        AddTH2("Cluster Size vs. Energy", "Energy [MeV]", "Cluster Size", Ebins, BinSettings(50), "clusterSize_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_clusterSize().at(i), f.TaggW());
        });

        AddTH1("Tagger Time - CB Average Time", "t [ns]", "#", TaggTime, "TaggTime",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
        });

        AddTH2("dEvEproton", "E [MeV]", "dE [MeV]", Ebins, vetoEbins, "dEvE",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.p_kinfitted().Energy() - ParticleTypeDatabase::Proton.Mass(), f.Tree.p_vetoE);
        });

    }

    void Fill(const Fill_t& f) const
    {
        h1.Fill(f);
        h2.Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        vector<TH1*> v;
        v.reserve(h1.size()+h2.size());
        for (auto& e : h1)
            v.emplace_back(e.h);
        for (auto& e: h2)
            v.emplace_back(e.h);
        return v;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts()
    {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb > 0.001", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .001; }},
                              {"KinFitProb > 0.02", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .02; }},
                              {"KinFitProb > 0.05", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .05; }}
                             });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"TreeFitProb > 0.001", [] (const Fill_t& f) { return f.Tree.treefit_probability > .001; }},
                              {"TreeFitProb > 0.01", [] (const Fill_t& f) { return f.Tree.treefit_probability > .01; }},
                              {"TreeFitProb > 0.02", [] (const Fill_t& f) { return f.Tree.treefit_probability > .02; }},
                              {"TreeFitProb > 0.05", [] (const Fill_t& f) { return f.Tree.treefit_probability > .05; }}
                          });

/*        auto antiPi0Cut = [] (const Fill_t& f, const double low = 102., const double high = 170.) {
            const interval<double> pion_cut(low, high);
            TLorentzVector pi0;
            const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};

            const auto sorted = std_ext::get_sorted_indices(f.Tree.photons_vetoE());

            for (const auto pi0_comb : pi0_combs) {
                pi0 = TLorentzVector(0., 0., 0., 0.);

                for (const auto idx : pi0_comb)
                    pi0 += TParticle(ParticleTypeDatabase::Photon, f.Tree.photons().at(sorted.at(idx)));

                // check anti pi^0 cut
                if (pion_cut.Contains(pi0.M()))
                    return false;
            }

            return true;
        };

        auto IM2d_lin_cut = [] (const Fill_t& f) {
            const double eeIM = im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted());

            return eeIM < (1.15*f.Tree.eta().M() - 170);
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"anti pi0", antiPi0Cut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"lin cut",  IM2d_lin_cut}
                          });
*/
//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"MM < 1030 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1030; }},
//                              {"MM < 1010 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1010; }},
//                              {"MM < 1000 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1000; }},
//                              {"MM < 990 MeV",  [] (const Fill_t& f) { return f.Tree.mm().M() < 990; }}
//                          });


        return cuts;
    }

};



struct MesonDalitzDecays_plot : Plotter {

    Hist_t::Tree_t Tree;

    using MCHist_t = MCTrue_Splitter<Hist_t>;
    cuttree::Tree_t<MCHist_t> hists;

    MesonDalitzDecays_plot(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {
        init_tree(input, Tree, "MesonDalitzDecays/tree");

        hists = cuttree::Make<MCHist_t>(HistFac);
    }

    static void init_tree(const WrapTFileInput& input, WrapTTree& tree, const string& name)
    {
        if (!input.GetObject(name, tree.Tree))
            throw Exception("Cannot find tree "+name);

        tree.LinkBranches();
    }

    virtual long long GetNumEntries() const override
    {
        return Tree.Tree->GetEntries();
    }

    virtual void ProcessEntry(const long long entry) override
    {
        Tree.Tree->GetEntry(entry);
        cuttree::Fill<MCHist_t>(hists, {Tree});
    }

};



AUTO_REGISTER_PLOTTER(MesonDalitzDecays_plot)
