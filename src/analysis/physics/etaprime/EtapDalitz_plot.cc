#include "etaprime_dalitz.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "base/Logger.h"
#include "base/interval.h"

#include "expconfig/ExpConfig.h"

#include "TCutG.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;


template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    const decltype(physics::EtapDalitz::makeChannels()) reaction_channels;

    constexpr static Color_t bkg_color = kGray+1;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
        reaction_channels(physics::EtapDalitz::makeChannels())
    {
        using histstyle::Mod_t;

        const Color_t sig_color = reaction_channels.channels.find(1)->second.color;
        const Color_t ref_color = reaction_channels.channels.find(2)->second.color;

        this->GetHist(0, "Data", Mod_t::MakeDataPoints(kBlack));
        this->GetHist(5, "D07", Mod_t::MakeDataPoints(kGray));
        this->GetHist(6, "D10", Mod_t::MakeDataPoints(kGray));
        this->GetHist(7, "D12", Mod_t::MakeDataPoints(kGray));

        this->GetHist(1, "Signal", Mod_t::MakeLine(sig_color, 2));
        this->GetHist(2, "Reference", Mod_t::MakeLine(ref_color, 2));

        // mctrue is never >= 4 (and < 9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(bkg_color, -1));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = unsigned(f.Tree.channel);

        using histstyle::Mod_t;

        auto get_bkg_name = [this] (const unsigned mctrue) {
            const auto entry = reaction_channels.channels.find(int(mctrue));

            if (entry != reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        auto get_color = [this] (const unsigned mctrue) -> short {
            const auto entry = reaction_channels.channels.find(int(mctrue));

            if (entry != reaction_channels.channels.end())
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
            this->GetHist(3).Fill(f);
            if (mctrue >= 10)
                this->GetHist(4).Fill(f);
        }

        // handle D07/D10/D12
        if (mctrue == 0)
            this->GetHist(4+f.Tree.beamtime).Fill(f);
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

template <typename T>
vector<size_t> get_sorted_indices(vector<T> vec)
{
    vector<size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
              [vec] (size_t i, size_t j) {
        return vec[i] > vec[j];
    });

    return p;
}

double im_ee(vector<double> vetoE, vector<TLorentzVector> photons)
{
    const auto leptons = get_sorted_indices(vetoE);

    return (photons.at(leptons[0]) + photons.at(leptons[1])).M();
}


template <typename Tree_t>
struct Hist_t {

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

    static constexpr double binScale = 1.;

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    HistogramFactory HistFac;

    Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf) {}

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

    cuttree::Cuts_t<Fill_t> GetCuts();

    struct TreeCuts {

        struct antiPi0Cut {
            const double low;
            const double high;
            antiPi0Cut(const double l = 102., const double h = 170.) : low(l), high(h) {}

            bool operator() (const Fill_t& f) const {
                const interval<double> pion_cut(low, high);
                TLorentzVector pi0;
                const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};

                const auto photons = f.Tree.photons();
                const auto sorted = get_sorted_indices(f.Tree.photons_vetoE());

                for (const auto pi0_comb : pi0_combs) {
                    pi0 = TLorentzVector(0., 0., 0., 0.);

                    for (const auto idx : pi0_comb)
                        pi0 += TParticle(ParticleTypeDatabase::Photon, photons.at(sorted.at(idx)));

                    // check anti pi^0 cut
                    if (pion_cut.Contains(pi0.M()))
                        return false;
                }

                return true;
            }
        };

        static bool distinctPIDCut(const Fill_t& f) noexcept {
            const auto channels = f.Tree.photons_vetoChannel();
            const auto idx = get_sorted_indices(f.Tree.photons_vetoE());

            return channels.at(idx[0]) != channels.at(idx[1]);
        }

        struct freeZ_vertexCut {
            const double high;
            freeZ_vertexCut(const double h = 6.) : high(h) {}

            bool operator() (const Fill_t& f) const {
                return f.Tree.kinfit_freeZ_ZVertex < high;
            }
        };

        struct treefit_vertexCut {
            const double low;
            const double high;
            treefit_vertexCut(const double l = -7., const double h = 7.) : low(l), high(h) {}

            bool operator() (const Fill_t& f) const {
                return f.Tree.treefit_ZVertex > low && f.Tree.treefit_ZVertex < high;
            }
        };

        static bool pid_cut(const Fill_t& f, const double threshold) {
            const auto vetos = f.Tree.photons_vetoE();
            const auto idx = get_sorted_indices(vetos);

            return vetos.at(idx[0]) > threshold && vetos.at(idx[1]) > threshold;
        }

        static bool allFS_CB(const Fill_t& f) noexcept {
            size_t nCB = 0;
            for (const auto& d : f.Tree.photons_detector())
                if (d == 1)
                    nCB++;

            if (nCB < 3)
                return false;
            return true;
        }
    };
};

// define the structs containing the histograms and the cuts
struct SigHist_t : Hist_t<physics::EtapDalitz::SigTree_t> {

    using Tree_t = physics::EtapDalitz::SigTree_t;
    using Fill_t = Hist_t<Tree_t>::Fill_t;

    const BinSettings Ebins    = Bins(1200, 0, 1200);

    const BinSettings Chi2Bins = BinSettings(250, 0, 25);
    const BinSettings probbins = BinSettings(250, 0,  1);

    const BinSettings IMbins   = Bins(1200,   0, 1200);
    const BinSettings MMbins   = Bins(1200, 400, 1600);

    const BinSettings TaggChBins = BinSettings(47);

    const BinSettings TaggTime   = BinSettings(240, -30, 30);
    const BinSettings CoplBins   = Bins(300, 0, 30);

    const BinSettings zVertex    = Bins(100, -15, 15);
    const BinSettings free_vz    = Bins(400, -40, 40);


    //HistogramFactory HistFac;

    SigHist_t(const HistogramFactory& hf, cuttree::TreeInfo_t treeInfo) : Hist_t(hf, treeInfo) {

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

        AddTH1("3 photon IM", "3#gamma IM [MeV]", "#", IMbins, "etapIM",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap().M(), f.TaggW());
        });

        AddTH1("3 photon IM kinfitted",  "3#gamma IM fit [MeV]", "#", IMbins, "etapIM_kinfitted",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap_kinfit().M(), f.TaggW());
        });

        AddTH1("Missing Mass", "MM [MeV]", "", MMbins, "mm",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm().M(), f.TaggW());
        });

        AddTH1("Z Vertex Kinfit", "z [cm]", "#", zVertex, "v_z_kinfit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex free Z Kinfit", "z [cm]", "#", free_vz, "v_z_kinfit_freeZ",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_freeZ_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex Treefit", "z [cm]", "#", zVertex, "v_z_treefit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex free Z Treefit", "z [cm]", "#", zVertex, "v_z_treefit_freeZ",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_freeZ_ZVertex, f.TaggW());
        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()), f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()), f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(f.Tree.photons_vetoE(), f.Tree.photons_kinfitted()));
//        });

//        AddTH2("Cluster Size vs. Energy", "Energy [MeV]", "Cluster Size", Ebins, BinSettings(50), "clusterSize_E",
//               [] (TH2D* h, const Fill_t& f) {
//            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
//                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_clusterSize().at(i), f.TaggW());
//        });

        AddTH1("Effective Cluster Radius", "R", "#", BinSettings(500, 0, 50), "clusterRadius",
               [] (TH1D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_effect_radius().at(i), f.TaggW());
        });

        AddTH1("Lateral Moment", "L", "#", BinSettings(200, 0, 1), "lateralMoment",
               [] (TH1D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons_lat_moment().at(i), f.TaggW());
        });

        AddTH2("Effective Cluster Radius vs. Energy", "Energy [MeV]", "R", Bins(300, 0, 1200), BinSettings(200, 0, 50), "clusterRadius_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_effect_radius().at(i), f.TaggW());
        });

        AddTH2("Lateral Moment vs. Energy", "Energy [MeV]", "L", Bins(300, 0, 1200), BinSettings(100, 0, 1), "lateralMoment_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons_lat_moment().at(i), f.TaggW());
        });

//        AddTH2("Lateral Moment vs. Effective Cluster Radius", "R", "L", BinSettings(500, 0, 50), BinSettings(200, 0, 1), "lateralMoment_clusterRadius",
//               [] (TH2D* h, const Fill_t& f) {
//            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
//                h->Fill(f.Tree.photons_effect_radius().at(i), f.Tree.photons_lat_moment().at(i), f.TaggW());
//        });

//        AddTH1("Tagger Time - CB Average Time", "t [ns]", "#", TaggTime, "TaggTime",
//               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
//        });


        AddTH1("TOF TAPS photon", "TOF [ns]", "#", TaggTime, "TOF_gTAPS",
               [] (TH1D* h, const Fill_t& f) {
            const auto idx = get_sorted_indices(f.Tree.photons_vetoE());
            if (f.Tree.photons_detector().at(idx[2]) != 2)
                return;
            h->Fill(f.Tree.photons_Time().at(idx[2]));
        });

    }

    static TCutG* makeEffectiveRadiusCut()
    {
        TCutG* c = new TCutG("EffectiveRadiusCut", 5);
        c->SetPoint(0, 1200.,  6.);
        c->SetPoint(1,  800.,  8.);
        c->SetPoint(2,  800., 13.);
        c->SetPoint(3, 1200., 13.);
        c->SetPoint(3, 1200.,  6.);
        return c;
    }

    static TCutG* makeLateralMomentCut()
    {
        TCutG* c = new TCutG("LateralMomentCut", 7);
        c->SetPoint(0,  80., 1.);
        c->SetPoint(1,  40.,  .9);
        c->SetPoint(2, 600.,  .35);
        c->SetPoint(3, 600., 0.);
        c->SetPoint(4,   0., 0.);
        c->SetPoint(5,   0., 1.);
        c->SetPoint(6,  80., 1.);
        return c;
    }

    static TCutG* makeSmallLateralMomentCut()
    {
        TCutG* c = new TCutG("SmallLateralMomentCut", 5);
        c->SetPoint(0, 100., 1.);
        c->SetPoint(1,  80.,  .85);
        c->SetPoint(2,   0.,  .85);
        c->SetPoint(3,   0., 1.);
        c->SetPoint(4, 100., 1.);
        return c;
    }

    static TCutG* effectiveRadiusCut;
    static TCutG* lateralMomentCut;
    static TCutG* smallLateralMomentCut;

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts()
    {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb > 0.001", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .001; }},
                              {"KinFitProb > 0.02", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .02; }},
                              {"KinFitProb > 0.05", [] (const Fill_t& f) { return f.Tree.kinfit_probability > .05; }}
                          });

        auto eff_radius_cut = [] (const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (effectiveRadiusCut->IsInside(f.Tree.photons().at(i).Energy(), f.Tree.photons_effect_radius().at(i)))
                    return false;
            return true;
        };

        auto lat_moment_cut = [] (const Fill_t& f, const TCutG* const cut) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (cut->IsInside(f.Tree.photons().at(i).Energy(), f.Tree.photons_lat_moment().at(i)))
                    return false;
            return true;
        };

        auto lateral_moment = [] (const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (f.Tree.photons_lat_moment().at(i) > .95)
                    return false;
            return true;
        };

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"distinct PID elements", TreeCuts::distinctPIDCut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"anti pi0", TreeCuts::antiPi0Cut()}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"free vz cut", TreeCuts::freeZ_vertexCut()},
                              {"treefit vz cut", TreeCuts::treefit_vertexCut()}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"effective radius", eff_radius_cut}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"lateral moment", [&lat_moment_cut] (const Fill_t& f) {
                                   return lat_moment_cut(f, lateralMomentCut);
                               }},
                              {"small lateral moment", [&lat_moment_cut] (const Fill_t& f) {
                                   return lat_moment_cut(f, smallLateralMomentCut);
                               }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"lateral moment < .95", lateral_moment}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"MM < 1030 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1030; }},
                              {"MM < 1010 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1010; }},
                              {"MM < 1000 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1000; }},
                              {"MM < 990 MeV",  [] (const Fill_t& f) { return f.Tree.mm().M() < 990; }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"PID e^{#pm} > .4 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .4); }},
                              {"PID e^{#pm} > .5 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .5); }},
                              {"PID e^{#pm} > .6 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .6); }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"allFS in CB", TreeCuts::allFS_CB}
                          });

        return cuts;
    }

};

struct RefHist_t : Hist_t<physics::EtapDalitz::RefTree_t> {

    using Tree_t = physics::EtapDalitz::RefTree_t;
    using Fill_t = Hist_t<Tree_t>::Fill_t;

    const BinSettings Chi2Bins = BinSettings(250, 0, 25);
    const BinSettings probbins = BinSettings(250, 0,  1);

    const BinSettings IMbins   = Bins(1200,   0, 1200);

    const BinSettings zVertex  = Bins(100, -15, 15);

    RefHist_t(const HistogramFactory& hf, cuttree::TreeInfo_t treeInfo) : Hist_t(hf, treeInfo) {

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

        AddTH1("2 photon IM", "2#gamma IM [MeV]", "#", IMbins, "etapIM",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap().M(), f.TaggW());
        });

        AddTH1("2 photon IM kinfitted",  "2#gamma IM fit [MeV]", "#", IMbins, "etapIM_kinfitted",
               [] (TH1D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap_kinfit().M(), f.TaggW());
        });

        AddTH1("Z Vertex Kinfit", "z [cm]", "#", zVertex, "v_z_kinfit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_ZVertex, f.TaggW());
        });

        AddTH1("Z Vertex Treefit", "z [cm]", "#", zVertex, "v_z_treefit",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_ZVertex, f.TaggW());
        });

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
                              {"treefit vz cut", TreeCuts::treefit_vertexCut()}
                          });

        return cuts;
    }

};



template <typename Hist_t>
struct EtapDalitz_plot : Plotter {

    typename Hist_t::Tree_t Tree;

    using MCHist_t = MCTrue_Splitter<Hist_t>;
    cuttree::Tree_t<MCHist_t> cuttree_hists;

    EtapDalitz_plot(const string& tag, const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {
        const auto tree_base = opts->Get<string>("Tree", "EtapDalitz");

        init_tree(input, Tree, tree_base + "/" + tag);
        cuttree_hists = cuttree::Make<MCHist_t>(HistFac);
    }

    static void init_tree(const WrapTFileInput& input, WrapTTree& tree, const string& name)
    {
        if (!input.GetObject(name, tree.Tree))
            throw Exception("Can't find tree " + name);

        tree.LinkBranches();
    }

    virtual long long GetNumEntries() const override
    {
        return Tree.Tree->GetEntries();
    }

    virtual void ProcessEntry(const long long entry) override
    {
        Tree.Tree->GetEntry(entry);
        cuttree::Fill<MCHist_t>(cuttree_hists, {Tree});
    }

};


struct EtapDalitz_plot_Sig : EtapDalitz_plot<SigHist_t> {

    EtapDalitz_plot_Sig(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapDalitz_plot("signal", name, input, opts)
    {
    }
};


struct EtapDalitz_plot_Ref : EtapDalitz_plot<RefHist_t> {

    EtapDalitz_plot_Ref(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapDalitz_plot("ref", name, input, opts)
    {
    }
};


TCutG* SigHist_t::effectiveRadiusCut = SigHist_t::makeEffectiveRadiusCut();
TCutG* SigHist_t::lateralMomentCut = SigHist_t::makeLateralMomentCut();
TCutG* SigHist_t::smallLateralMomentCut = SigHist_t::makeSmallLateralMomentCut();


AUTO_REGISTER_PLOTTER(EtapDalitz_plot_Sig)
AUTO_REGISTER_PLOTTER(EtapDalitz_plot_Ref)
