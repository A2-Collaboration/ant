#include "etaprime_dalitz.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "base/Logger.h"
#include "base/interval.h"
#include "base/std_ext/vector.h"

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
        this->GetHist(5, "D07", Mod_t::MakeDataPoints(kGreen-8));
        this->GetHist(6, "D10", Mod_t::MakeDataPoints(kOrange-8));
        this->GetHist(7, "D12", Mod_t::MakeDataPoints(kBlue-8));

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
            const auto entry = reaction_channels.channels.find(mctrue);

            if (entry != reaction_channels.channels.end())
                return entry->second.name;

            return string("Unknown Decay");
        };

        auto get_color = [this] (const unsigned mctrue) -> short {
            const auto entry = reaction_channels.channels.find(mctrue);

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

vector<double> get_veto_energies(vector<TSimpleParticle> particles)
{
    vector<double> veto_energies;
    for (const auto& p : particles)
        veto_energies.emplace_back(p.VetoE);

    return veto_energies;
}

vector<size_t> get_sorted_indices_vetoE(vector<TSimpleParticle> particles)
{
    return std_ext::get_sorted_indices_desc(get_veto_energies(particles));
}

double im_ee(vector<TSimpleParticle> photons)
{
    const auto leptons = get_sorted_indices_vetoE(photons);

    return (photons.at(leptons[0]) + photons.at(leptons[1])).M();
}

double im_ee(vector<double> vetoE, vector<TLorentzVector> photons)
{
    const auto leptons = std_ext::get_sorted_indices_desc(vetoE);

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
    const bool isLeaf;

    static BinSettings Bins(const unsigned bins, const double min, const double max) {
        return BinSettings(unsigned(bins*binScale), min, max);
    }

    HistMgr<TH1D> h1;
    HistMgr<TH2D> h2;

    HistogramFactory HistFac;

    Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t treeInfo) :
        isLeaf(treeInfo.nDaughters == 0),
        HistFac(hf)
    {}

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

        static constexpr bool prob_cut(const double prob_val, const double prob_thresh) {
            return !std::isnan(prob_val) && prob_val > prob_thresh;
        }

        // in a few very rare cases the probability is 1.0 and the chi2 is nan
        // the fit most likely failed but responded to be successful, take this into account
        static constexpr bool prob_cut(const double prob_val, const double prob_thresh, const double chi2) {
            return !std::isnan(prob_val) && prob_val > prob_thresh && !std::isnan(chi2);
        }

        struct antiPi0Cut {
            const double low;
            const double high;
            antiPi0Cut(const double l = 102., const double h = 170.) : low(l), high(h) {}

            bool operator() (const Fill_t& f) const {
                const interval<double> pion_cut(low, high);
                TLorentzVector pi0;
                const std::vector<std::array<size_t, 2>> pi0_combs = {{0, 2}, {1, 2}};

                const auto photons = f.Tree.photons();
                const auto sorted = get_sorted_indices_vetoE(photons);

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

        static bool antiPi0Cut_sigma(const Fill_t& f, const double Nsigma) {
            // values obtained from kinematically unfitted pi0 peak from subIM combs of eta'
            static constexpr double mean = 136.3;
            static constexpr double sigma = 12.01;

            const double cut_sigma = Nsigma*sigma;
            const auto cut = antiPi0Cut(mean-cut_sigma, mean+cut_sigma);
            return cut(f);
        }

        static bool distinctPIDCut(const Fill_t& f) noexcept {
            const auto channels = f.Tree.photons_vetoChannel();
            const auto idx = get_sorted_indices_vetoE(f.Tree.photons());

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
            const auto vetos = get_veto_energies(f.Tree.photons());
            const auto idx = std_ext::get_sorted_indices_desc(vetos);

            return vetos.at(idx[0]) > threshold && vetos.at(idx[1]) > threshold;
        }

        static bool allFS_CB(const Fill_t& f) noexcept {
            size_t nCB = count_if(f.Tree.photons_detector().begin(), f.Tree.photons_detector().end(),
                                  [](const int d){ return d == 1; });

            if (nCB < 3)
                return false;
            return true;
        }

        static bool lateral_moment(const Fill_t& f, const double threshold) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (f.Tree.photons_lat_moment().at(i) > threshold)
                    return false;
            return true;
        }

        static bool eff_radius_2d_cut(const Fill_t& f, const TCutG* const cut) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (cut->IsInside(f.Tree.photons().at(i).Energy(), f.Tree.photons_effect_radius().at(i)))
                    return false;
            return true;
        }

        static bool lat_moment_2d_cut(const Fill_t& f, const TCutG* const cut) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (cut->IsInside(f.Tree.photons().at(i).Energy(), f.Tree.photons_lat_moment().at(i)))
                    return false;
            return true;
        }

        static bool cluster_size_2d_cut(const Fill_t& f, const TCutG* const cut) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                if (cut->IsInside(f.Tree.photons().at(i).Energy(), f.Tree.photons().at(i).ClusterSize))
                    return false;
            return true;
        }

        static bool discarded_energy(const Fill_t& f, const double threshold) {
            return f.Tree.DiscardedEk <= threshold;
        }

        static bool hard_select(const Fill_t& f) noexcept {
            return allFS_CB(f) && distinctPIDCut(f) && discarded_energy(f, 0.);
        }

        static bool im900(const Fill_t& f) noexcept {
            return f.Tree.etap_kinfit().M() > 900.;
        }

        static constexpr bool do_nothing(const Fill_t&) noexcept {
            return true;
        }
    };

    // results from a PCA with the variables for lateral moment, cluster size, and cluster energy
    struct PCA_ClusterShape_t {
        //
        // Static data variables
        //
        static constexpr int gNVariables = 3;

        // Assignment of eigenvector matrix.
        // Elements are stored row-wise, that is
        //    M[i][j] = e[i * nVariables + j]
        // where i and j are zero-based
        static constexpr double gEigenVectors[] = {
            0.401069,
            0.807083,
            0.433315,
            0.716885,
            0.0179446,
            -0.69696,
            0.57028,
            -0.590166,
            0.571389};

        // Assignment to eigen value vector. Zero-based.
        static constexpr double gEigenValues[] = {
            0.566304,
            0.346853,
            0.0868428
        };

        // Assignment to mean value vector. Zero-based.
        static constexpr double gMeanValues[] = {
            0.819355,
            8.7192,
            347.709
        };

        // Assignment to sigma value vector. Zero-based.
        static constexpr double gSigmaValues[] = {
            0.191971,
            3.84283,
            208.744
        };

        //
        // The function   void X2P(Double_t *x, Double_t *p)
        //
        static void X2P(Double_t *x, Double_t *p) {
            for (Int_t i = 0; i < gNVariables; i++) {
                p[i] = 0;
                for (Int_t j = 0; j < gNVariables; j++)
                    p[i] += (x[j] - gMeanValues[j])
                            * gEigenVectors[j *  gNVariables + i] / gSigmaValues[j];

            }
        }

        //
        // The function   void P2X(Double_t *p, Double_t *x, Int_t nTest)
        //
        static void P2X(Double_t *p, Double_t *x, Int_t nTest) {
            for (Int_t i = 0; i < gNVariables; i++) {
                x[i] = gMeanValues[i];
                for (Int_t j = 0; j < nTest; j++)
                    x[i] += p[j] * gSigmaValues[i]
                            * gEigenVectors[i *  gNVariables + j];
            }
        }

    };
};

// make the linker happy
template <typename Tree_t>
constexpr int Hist_t<Tree_t>::PCA_ClusterShape_t::gNVariables;
template <typename Tree_t>
constexpr double Hist_t<Tree_t>::PCA_ClusterShape_t::gEigenVectors[];
template <typename Tree_t>
constexpr double Hist_t<Tree_t>::PCA_ClusterShape_t::gEigenValues[];
template <typename Tree_t>
constexpr double Hist_t<Tree_t>::PCA_ClusterShape_t::gMeanValues[];
template <typename Tree_t>
constexpr double Hist_t<Tree_t>::PCA_ClusterShape_t::gSigmaValues[];


template <typename Tree_t>
struct q2Hist_t {

    // use needed types/structs which are defined in Hist_t
    using Fill_t = typename Hist_t<Tree_t>::Fill_t;
    template <typename Hist>
    using HistFiller_t = typename Hist_t<Tree_t>::template HistFiller_t<Hist>;

    struct q2_params_t {
        static constexpr double bin_width = 50.;
        static constexpr double max_value = 900.;
    };

    template <typename Hist>
    struct q2HistMgr : std::list<HistFiller_t<Hist>> {

        using list<HistFiller_t<Hist>>::list;

        void Fill(const Fill_t& data) const
        {
            const auto imee = im_ee(get_veto_energies(data.Tree.photons()), data.Tree.photons_kinfitted());

            // make sure the momentum transfer has physical reasonable values
            if (imee < 0. || !isfinite(imee))
                return;

            if (imee > q2_params_t::max_value)
                return;

            size_t idx = imee/q2_params_t::bin_width;
            auto it = this->begin();
            std::advance(it, idx);
            it->Fill(data);
        }
    };

    q2HistMgr<TH1D> q2_hists;

    HistogramFactory q2_hf;

    q2Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t) :
        // don't create subfolder, doesn't work with hstack, seems too complicated to modify it to get it work atm...
        //q2_hf(HistogramFactory("q2_bins", hf))
        q2_hf(hf)
    {
        LOG_N_TIMES(1, INFO) << "Use fixed binning for q2 hists with bin width = " << q2_params_t::bin_width;
        for (double q2 = 0.; q2 < q2_params_t::max_value; q2 += q2_params_t::bin_width) {
            stringstream ss_title;
            stringstream ss_id;
            ss_title << "IMee " << q2 << " to " << q2+q2_params_t::bin_width << " MeV";
            ss_id << "imee_" << q2 << "_" << q2+q2_params_t::bin_width;
            q2_hists.emplace_back(HistFiller_t<TH1D>(
                                      q2_hf.makeTH1D(ss_title.str(), "IM(eeg) [MeV]", "#", BinSettings(240, 0, 1200), ss_id.str()),
                                      [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.etap_kinfit().M(), f.TaggW()); }));
        }
    }

    void Fill(const Fill_t& f) const
    {
        q2_hists.Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        vector<TH1*> v;
        for (auto& e: q2_hists)
            v.emplace_back(e.h);
        return v;
    }

    cuttree::Cuts_t<Fill_t> GetCuts();
};

// variable bin width variant of the q2Hist_t struct defined above
template <typename Tree_t>
struct q2Hist_var_t {

    // use needed types/structs which are defined in Hist_t
    using Fill_t = typename Hist_t<Tree_t>::Fill_t;
    template <typename Hist>
    using HistFiller_t = typename Hist_t<Tree_t>::template HistFiller_t<Hist>;

    struct q2_params_t {
        static constexpr double min_value = 0.;
        static constexpr double max_value = 900.;
        static vector<double> bin_widths;
    };

    template <typename Hist>
    struct q2HistMgr : std::list<HistFiller_t<Hist>> {

        using list<HistFiller_t<Hist>>::list;

        void Fill(const Fill_t& data) const
        {
            const double imee = im_ee(get_veto_energies(data.Tree.photons()), data.Tree.photons_kinfitted());

            // make sure the momentum transfer has physical reasonable values
            if (imee < 0. || !isfinite(imee))
                return;

            if (imee > q2_params_t::max_value)
                return;

            auto it = this->begin();
            auto bins = q2_params_t::bin_widths.begin();
            // advance histograms iterator as long as the accumulated q2 value of the bin widths is below the current value
            for (double q2 = 0.; q2 < imee; q2 += *bins++, ++it);
            it->Fill(data);
        }
    };

    q2HistMgr<TH1D> q2_hists;

    HistogramFactory q2_hf;

    q2Hist_var_t(const HistogramFactory& hf, cuttree::TreeInfo_t, const vector<double> bins) :
        // don't create subfolder, doesn't work with hstack, seems too complicated to modify it to get it work atm...
        //q2_hf(HistogramFactory("q2_bins", hf))
        q2_hf(hf)
    {
        LOG_N_TIMES(1, INFO) << "Use variable q2 hists";
        q2_params_t::bin_widths = bins;
        double bin_start = q2_params_t::min_value;
        auto it = q2_params_t::bin_widths.begin();

        while (bin_start < q2_params_t::max_value) {
            // sanity check to make sure enough bin widths are provided to cover the whole region
            if (it == q2_params_t::bin_widths.end()) {
                LOG(ERROR) << "Not enough bins provided, max q^2 value of " << q2_params_t::max_value << " not reached";
                LOG(INFO) << "Given bin widths only cover the region up until " << bin_start;;
                throw runtime_error("Not enough bins provided");
            }

            // create the histograms for the bins
            double q2 = bin_start + *it++;
            stringstream ss_title;
            stringstream ss_id;
            ss_title << "IMee " << bin_start << " to " << q2 << " MeV";
            ss_id << "imee_" << bin_start << "_" << q2;
            q2_hists.emplace_back(HistFiller_t<TH1D>(
                                      q2_hf.makeTH1D(ss_title.str(), "IM(eeg) [MeV]", "#", BinSettings(240, 0, 1200), ss_id.str()),
                                      [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.etap_kinfit().M(), f.TaggW()); }));

            LOG_N_TIMES(bins.size(), INFO) << "Created q2 hist for range " << ss_title.str();
            bin_start = q2;
        }
    }

    void Fill(const Fill_t& f) const
    {
        q2_hists.Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        vector<TH1*> v;
        for (auto& e: q2_hists)
            v.emplace_back(e.h);
        return v;
    }

    cuttree::Cuts_t<Fill_t> GetCuts();
};

// help the linker figuring out the reference to the bin widths vector defined and used above
template<typename Tree_t>
vector<double> q2Hist_var_t<Tree_t>::q2_params_t::bin_widths;

// define the structs containing the histograms and the cuts
struct SigHist_t : Hist_t<physics::EtapDalitz::SigTree_t>, q2Hist_var_t<physics::EtapDalitz::SigTree_t> {

    using Tree_t = physics::EtapDalitz::SigTree_t;
    using Fill_t = Hist_t<Tree_t>::Fill_t;

    void Fill(const Fill_t& f) const
    {
        Hist_t::Fill(f);
        q2Hist_var_t::Fill(f);
    }

    std::vector<TH1*> GetHists() const
    {
        auto v = Hist_t::GetHists();
        std_ext::concatenate(v, q2Hist_var_t::GetHists());
        return v;
    }

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

    SigHist_t(const HistogramFactory& hf, cuttree::TreeInfo_t treeInfo) :
        Hist_t(hf, treeInfo),
        q2Hist_var_t(hf, treeInfo, {80., 60., 50., 50., 80., 100., 80., 60., 50., 50., 50., 50., 50., 50., 50.})
    {

/*
        AddTH1("PCA_1 shower shape", "pca1", "#", Bins(1000,-10,10), "pca1ShowerShape",
               [] (TH1D* h, const Fill_t& f) {
            vector<double> x;
            vector<double> p(PCA_ClusterShape_t::gNVariables);
            for (unsigned i = 0; i < f.Tree.photons().size(); i++) {
                x.emplace_back(f.Tree.photons_lat_moment().at(i));
                x.emplace_back(f.Tree.photons().at(i).ClusterSize);
                x.emplace_back(f.Tree.photons().at(i).Energy());
                PCA_ClusterShape_t::X2P(&x[0], &p[0]);
                x.clear();
                h->Fill(p[0], f.TaggW());
            }
        });
        AddTH1("PCA_2 shower shape", "pca2", "#", Bins(1000,-10,10), "pca2ShowerShape",
               [] (TH1D* h, const Fill_t& f) {
            vector<double> x;
            vector<double> p(PCA_ClusterShape_t::gNVariables);
            for (unsigned i = 0; i < f.Tree.photons().size(); i++) {
                x.emplace_back(f.Tree.photons_lat_moment().at(i));
                x.emplace_back(f.Tree.photons().at(i).ClusterSize);
                x.emplace_back(f.Tree.photons().at(i).Energy());
                PCA_ClusterShape_t::X2P(&x[0], &p[0]);
                x.clear();
                h->Fill(p[1], f.TaggW());
            }
        });
        AddTH1("PCA_3 shower shape", "pca3", "#", Bins(1000,-10,10), "pca3ShowerShape",
               [] (TH1D* h, const Fill_t& f) {
            vector<double> x;
            vector<double> p(PCA_ClusterShape_t::gNVariables);
            for (unsigned i = 0; i < f.Tree.photons().size(); i++) {
                x.emplace_back(f.Tree.photons_lat_moment().at(i));
                x.emplace_back(f.Tree.photons().at(i).ClusterSize);
                x.emplace_back(f.Tree.photons().at(i).Energy());
                PCA_ClusterShape_t::X2P(&x[0], &p[0]);
                x.clear();
                h->Fill(p[2], f.TaggW());
            }
        });
        AddTH2("Shower Shape: PCA2 vs. PCA1", "pca1", "pca2", Bins(500,-10,10), Bins(500,-10,10), "pca1pca2",
               [] (TH2D* h, const Fill_t& f) {
            vector<double> x;
            vector<double> p(PCA_ClusterShape_t::gNVariables);
            for (unsigned i = 0; i < f.Tree.photons().size(); i++) {
                x.emplace_back(f.Tree.photons_lat_moment().at(i));
                x.emplace_back(f.Tree.photons().at(i).ClusterSize);
                x.emplace_back(f.Tree.photons().at(i).Energy());
                PCA_ClusterShape_t::X2P(&x[0], &p[0]);
                x.clear();
                h->Fill(p[0], p[1], f.TaggW());
            }
        });
        return;

        AddTH1("KinFitChi2", "#chi^{2}", "#", Chi2Bins, "KinFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.kinfit_chi2, f.TaggW());
        });

        AddTH1("TreeFitChi2", "#chi^{2}", "#", Chi2Bins, "TreeFitChi2",
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.treefit_chi2, f.TaggW());
        });
*/
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
               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.mm, f.TaggW());
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

//        AddTH1("Discarded Energy", "discarded Ek [MeV]", "#", Bins(500, 0, 100), "discardedEk",
//               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.DiscardedEk, f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g)", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons()), f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap().M(), im_ee(f.Tree.photons()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit",
//               [] (TH2D* h, const Fill_t& f) {
//            h->Fill(f.Tree.etap_kinfit().M(), im_ee(get_veto_energies(f.Tree.photons()), f.Tree.photons_kinfitted()), f.TaggW());
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit prompt", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_prompt",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() > 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(get_veto_energies(f.Tree.photons()), f.Tree.photons_kinfitted()));
//        });

//        AddTH2("IM(e+e-) vs. IM(e+e-g) fit random", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(600, 0, 1200), BinSettings(500, 0, 1000), "IM2d_fit_random",
//               [] (TH2D* h, const Fill_t& f) {
//            if (f.TaggW() < 0)
//                h->Fill(f.Tree.etap_kinfit().M(), im_ee(get_veto_energies(f.Tree.photons()), f.Tree.photons_kinfitted()));
//        });

        AddTH2("Cluster Size vs. Energy", "Energy [MeV]", "Cluster Size", Ebins, BinSettings(50), "clusterSize_E",
               [] (TH2D* h, const Fill_t& f) {
            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
                h->Fill(f.Tree.photons().at(i).Energy(), f.Tree.photons().at(i).ClusterSize, f.TaggW());
        });

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

//        AddTH2("Lateral Moment vs. Effective Cluster Radius", "R", "L", BinSettings(500, 0, 50), BinSettings(200, 0, 1), "lateralMoment_clusterRadius",
//               [] (TH2D* h, const Fill_t& f) {
//            for (unsigned i = 0; i < f.Tree.photons().size(); i++)
//                h->Fill(f.Tree.photons_effect_radius().at(i), f.Tree.photons_lat_moment().at(i), f.TaggW());
//        });

//        AddTH1("Tagger Time - CB Average Time", "t [ns]", "#", TaggTime, "TaggTime",
//               [] (TH1D* h, const Fill_t& f) { h->Fill(f.Tree.TaggT - f.Tree.CBAvgTime);
//        });


//        AddTH1("TOF TAPS photon", "TOF [ns]", "#", TaggTime, "TOF_gTAPS",
//               [] (TH1D* h, const Fill_t& f) {
//            const auto idx = get_sorted_indices_vetoE(f.Tree.photons());
//            if (f.Tree.photons_detector().at(idx[2]) != 2)
//                return;
//            h->Fill(f.Tree.photons().at(idx[2]).Time);
//        });

//        if (!isLeaf)
//            return;

        AddTH2("IM(e+e-) vs. IM(e+e-g) [TFF]", "IM(e+e-g) [MeV]", "IM(e+e-) [MeV]", BinSettings(240, 0, 1200), BinSettings(100, 0, 1000), "TFFextract",
               [] (TH2D* h, const Fill_t& f) {
            h->Fill(f.Tree.etap_kinfit().M(), im_ee(get_veto_energies(f.Tree.photons()), f.Tree.photons_kinfitted()), f.TaggW());
        });

        if (isLeaf)
            return;

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

    }

//    static TCutG* makeEffectiveRadiusCut()
//    {
//        TCutG* c = new TCutG("EffectiveRadiusCut", 5);
//        c->SetPoint(0, 1200.,  6.);
//        c->SetPoint(1,  800.,  8.);
//        c->SetPoint(2,  800., 13.);
//        c->SetPoint(3, 1200., 13.);
//        c->SetPoint(4, 1200.,  6.);
//        return c;
//    }

    static TCutG* makeEffectiveRadiusCut()
    {
        TCutG* c = new TCutG("EffectiveRadiusCut", 5);
        c->SetPoint(0, 12., 15.);
        c->SetPoint(1, 50.,  8.);
        c->SetPoint(2, 40.,  5.);
        c->SetPoint(3, 12.,  4.);
        c->SetPoint(4, 12., 15.);
        return c;
    }

    static TCutG* makeBigEffectiveRadiusCut()
    {
        TCutG* c = new TCutG("BigEffectiveRadiusCut", 7);
        c->SetPoint(0, 12., 25.);
        c->SetPoint(1, 80., 15.);
        c->SetPoint(2, 70., 10.);
        c->SetPoint(3, 50.,  5.);
        c->SetPoint(4, 40.,  2.);
        c->SetPoint(5, 12.,  2.);
        c->SetPoint(6, 12., 25.);
        return c;
    }

    static TCutG* makeLatMomentCut()
    {
        TCutG* c = new TCutG("LatMomentCut_Rtest", 4);
        c->SetPoint(0, 60., 1.);
        c->SetPoint(1, 12.,  .9);
        c->SetPoint(2, 12., 1.);
        c->SetPoint(3, 60., 1.);
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

    static TCutG* makeClusterSizeCut()
    {
        TCutG* c = new TCutG("ClusterSizeCut", 9);
        c->SetPoint(0, 100., 0);
        c->SetPoint(1, 100., 2);
        c->SetPoint(2, 150., 4);
        c->SetPoint(3, 200., 5);
        c->SetPoint(4, 250., 6);
        c->SetPoint(5, 350., 8);
        c->SetPoint(6, 400., 9);
        c->SetPoint(7, 400., 0);
        c->SetPoint(0, 100., 0);
        return c;
    }

    static TCutG* makeTightClusterSizeCut()
    {
        TCutG* c = new TCutG("TightClusterSizeCut", 9);
        c->SetPoint(0, 100., 0);
        c->SetPoint(1, 100., 3);
        c->SetPoint(2, 150., 5);
        c->SetPoint(3, 200., 6);
        c->SetPoint(4, 250., 7);
        c->SetPoint(5, 350., 9);
        c->SetPoint(6, 400., 9);
        c->SetPoint(7, 400., 0);
        c->SetPoint(0, 100., 0);
        return c;
    }

    static TCutG* effectiveRadiusCut;
    static TCutG* bigEffectiveRadiusCut;
    static TCutG* latMomentCut;
    static TCutG* lateralMomentCut;
    static TCutG* smallLateralMomentCut;
    static TCutG* clusterSizeCut;
    static TCutG* tightClusterSizeCut;

    // Sig and Ref channel share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts()
    {

        using cuttree::MultiCut_t;

        cuttree::Cuts_t<Fill_t> cuts;

//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"allFS in CB", TreeCuts::allFS_CB}
//                          });
//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"distinct PID elements", TreeCuts::distinctPIDCut}
//                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"selection", TreeCuts::hard_select}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              //{"KinFitProb > 0.001", [] (const Fill_t& f) {
                              //     return TreeCuts::prob_cut(f.Tree.kinfit_probability, .001); }},
                              {"KinFitProb > 0.02", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .02,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.05", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .05,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.1", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .1,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.2", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .2,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.3", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .3,
                                   f.Tree.kinfit_chi2); }},
                              {"IM(e+e-g) > 900 MeV", TreeCuts::im900}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"cluster size", [] (const Fill_t& f) {
                                   return TreeCuts::cluster_size_2d_cut(f, clusterSizeCut);
                               }},
                              {"tight cluster size", [] (const Fill_t& f) {
                                  return TreeCuts::cluster_size_2d_cut(f, tightClusterSizeCut);
                              }},
                              {"!cluster size", [] (const Fill_t& f) {
                                   return !(TreeCuts::cluster_size_2d_cut(f, clusterSizeCut));
                               }},
                              {"!tight cluster size", [] (const Fill_t& f) {
                                  return !(TreeCuts::cluster_size_2d_cut(f, tightClusterSizeCut));
                              }},
                              {"nothing", TreeCuts::do_nothing},
                              {"IM(e+e-g) > 900 MeV", TreeCuts::im900}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"free vz cut 6", TreeCuts::freeZ_vertexCut(6)},
                              {"free vz cut 4", TreeCuts::freeZ_vertexCut(4)},
                              {"free vz cut 3", TreeCuts::freeZ_vertexCut(3)},
                              {"free vz cut 2", TreeCuts::freeZ_vertexCut(2)},
                              {"free vz cut 1", TreeCuts::freeZ_vertexCut(1)},
                              {"free vz cut 0", TreeCuts::freeZ_vertexCut(0)},
                              {"treefit vz cut", TreeCuts::treefit_vertexCut()}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"lateral moment < .98", [] (const Fill_t& f) { return TreeCuts::lateral_moment(f, .98); }},
                              {"lateral moment < .97", [] (const Fill_t& f) { return TreeCuts::lateral_moment(f, .97); }},
                              {"treefit vz cut tighter", TreeCuts::treefit_vertexCut(-6, 6)},
                              {"IM(e+e-g) > 900 MeV", TreeCuts::im900}
                          });
/*
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"discarded Ek <= 0 MeV",  [] (const Fill_t& f) { return TreeCuts::discarded_energy(f,  0.); }},
                              {"discarded Ek <= 20 MeV", [] (const Fill_t& f) { return TreeCuts::discarded_energy(f, 20.); }},
                              {"discarded Ek <= 40 MeV", [] (const Fill_t& f) { return TreeCuts::discarded_energy(f, 40.); }},
                              {"discarded Ek <= 60 MeV", [] (const Fill_t& f) { return TreeCuts::discarded_energy(f, 60.); }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              //{"anti pi0", TreeCuts::antiPi0Cut()}
                              {"anti #pi^{0} 2#sigma", [] (const Fill_t& f) { return TreeCuts::antiPi0Cut_sigma(f, 2); }},
                              {"anti #pi^{0} 3#sigma", [] (const Fill_t& f) { return TreeCuts::antiPi0Cut_sigma(f, 3); }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"free vz cut", TreeCuts::freeZ_vertexCut()},
                              {"treefit vz cut", TreeCuts::treefit_vertexCut()}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              //{"effective radius", eff_radius_cut}
                              {"effective radius", [] (const Fill_t& f) {
                                   return TreeCuts::eff_radius_2d_cut(f, effectiveRadiusCut);
                               }},
                              {"big effective radius", [] (const Fill_t& f) {
                                   return TreeCuts::eff_radius_2d_cut(f, bigEffectiveRadiusCut);
                               }},
                              {"lateral moment (Rtest)", [] (const Fill_t& f) {
                                  return TreeCuts::lat_moment_2d_cut(f, latMomentCut);
                              }}
                          });

//        cuts.emplace_back(MultiCut_t<Fill_t>{
//                              {"lateral moment", [] (const Fill_t& f) {
//                                   return TreeCuts::lat_moment_2d_cut(f, lateralMomentCut);
//                               }},
//                              {"small lateral moment", [] (const Fill_t& f) {
//                                   return TreeCuts::lat_moment_2d_cut(f, smallLateralMomentCut);
//                               }}
//                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              //{"lateral moment < .97", [] (const Fill_t& f) { return TreeCuts::lateral_moment(f, .97); }},
                              {"lateral moment < .98", [] (const Fill_t& f) { return TreeCuts::lateral_moment(f, .98); }},
                              //{"lateral moment < .99", [] (const Fill_t& f) { return TreeCuts::lateral_moment(f, .99); }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"MM < 1030 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1030; }},
                              //{"MM < 1010 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1010; }},
                              {"MM < 1000 MeV", [] (const Fill_t& f) { return f.Tree.mm().M() < 1000; }},
                              {"MM < 990 MeV",  [] (const Fill_t& f) { return f.Tree.mm().M() < 990; }}
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              //{"PID e^{#pm} > .4 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .4); }},
                              {"PID e^{#pm} > .5 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .5); }},
                              {"PID e^{#pm} > .6 MeV", [] (const Fill_t& f) { return TreeCuts::pid_cut(f, .6); }}
                          });
*/
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
                              {"KinFitProb > 0.001", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .001,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.02", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .02,
                                   f.Tree.kinfit_chi2); }},
                              {"KinFitProb > 0.05", [] (const Fill_t& f) {
                                   return TreeCuts::prob_cut(f.Tree.kinfit_probability, .05,
                                   f.Tree.kinfit_chi2); }}
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
TCutG* SigHist_t::bigEffectiveRadiusCut = SigHist_t::makeBigEffectiveRadiusCut();
TCutG* SigHist_t::latMomentCut = SigHist_t::makeLatMomentCut();
TCutG* SigHist_t::lateralMomentCut = SigHist_t::makeLateralMomentCut();
TCutG* SigHist_t::smallLateralMomentCut = SigHist_t::makeSmallLateralMomentCut();
TCutG* SigHist_t::clusterSizeCut = SigHist_t::makeClusterSizeCut();
TCutG* SigHist_t::tightClusterSizeCut = SigHist_t::makeTightClusterSizeCut();


AUTO_REGISTER_PLOTTER(EtapDalitz_plot_Sig)
AUTO_REGISTER_PLOTTER(EtapDalitz_plot_Ref)
