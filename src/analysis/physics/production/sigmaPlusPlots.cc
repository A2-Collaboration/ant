#include "sigmaPlus.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"
#include "base/std_ext/string.h"

#include "base/Logger.h"

#include "TH1D.h"

#include "plot/CutTree.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::plot;

using namespace std;

using WrapTree = sigmaPlus::PionProdTree;

OptionsPtr global_opts = nullptr;


class sigmaPlus_Plot: public Plotter {

protected:

    template<typename Hist_t>
    struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

        // Hist_t should have that type defined
        using Fill_t = typename Hist_t::Fill_t;



        MCTrue_Splitter(const HistogramFactory& histFac,
                        const cuttree::TreeInfo_t& treeInfo) :
            cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)

        {
            using histstyle::Mod_t;


            this->GetHist(sigmaPlus::settings_t::Index_Data,
                          data_name, Mod_t::MakeDataPoints(kBlack));
            this->GetHist(sigmaPlus::settings_t::Index_Signal,
                          "Sig",  Mod_t::MakeLine(kRed, 2));
            this->GetHist(sigmaPlus::settings_t::Index_MainBkg,
                          "MainBkg",  Mod_t::MakeLine(kGreen, 2));


            this->GetHist(sigmaPlus::settings_t::Index_SumMC,
                          "Sum_MC", Mod_t::MakeLine(kBlack, 1));
            this->GetHist(sigmaPlus::settings_t::Index_BkgMC,
                          "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

            this->GetHist(sigmaPlus::settings_t::Index_brokenTree,
                          "brokenTree", Mod_t::MakeLine(kGray,1,-1));
            this->GetHist(sigmaPlus::settings_t::Index_unregTree,
                          "untaggedTree", Mod_t::MakeLine(kGray,1,-1));



            for (auto i_offset = 0u ; i_offset < sigmaPlus::otherBackgrounds.size() ; ++i_offset)
            {
                const auto mctrue = sigmaPlus::settings_t::Index_Offset + i_offset;
                this->GetHist(mctrue,
                              physics::sigmaPlus::getOtherChannelNames(mctrue),
                              Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-10), 1, kGray+1) );
            }

        }

        void Fill(const Fill_t& f) {

            const unsigned mctrue = f.Tree.MCTrue;
            const Hist_t& hist = this->GetHist(mctrue);
            hist.Fill(f);

            // handle MC_all and MC_bkg
            if(mctrue>0) {
                this->GetHist(3).Fill(f);
                if(mctrue >= 8 || mctrue == 2)
                    this->GetHist(4).Fill(f);
            }
        }
    };

    struct SigmaK0Hist_t {

        using Tree_t  = physics::sigmaPlus::PionProdTree;
        using RTree_t = sigmaPlus::RecTree;


        struct Fill_t {
            const Tree_t&  Tree;
            const RTree_t& RTree;
            Fill_t(const Tree_t& t, const RTree_t& rtree) : Tree(t), RTree(rtree) {}

            double TaggW() const {
                return Tree.Tagg_W;
            }

            vector<TLorentzVector> get2G(const vector<TSimpleParticle>& photons) const
            {
                vector<TLorentzVector> acc;
                const auto& permutation(Tree.SIG_combination());
                for (size_t i = 0 ; i < permutation.size(); i+=2)
                {
                    TLorentzVector gg(photons.at(permutation.at(i)));
                    gg += photons.at(permutation.at(i+1));
                    acc.emplace_back(gg);
                }
                return acc;
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

        const BinSettings taggerBins = sigmaPlus_Plot::taggerBins;


        const BinSettings Ebins    = Bins(1000, 0, 1000);

        const BinSettings Chi2Bins = Bins(250, 0,   25);
        const BinSettings probbins = Bins(250, 0,   1);

        const BinSettings IMbins       = Bins(1000,  200, 1100);
        const BinSettings IMProtonBins = Bins(1000,  600, 1200);
        const BinSettings IM2g         = Bins(1000,    0,  360);

        const BinSettings pThetaBins = Bins( 200,  0,   80);
        const BinSettings pEbins     = Bins( 350,  0, 1200);

        HistogramFactory HistFac;

        void AddTH1(const string &title, const string &xlabel, const string &ylabel,
                    const BinSettings &bins, const string &name,
                    const bool sumw2, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name, sumw2),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel,
                    const BinSettings &xbins, const BinSettings& ybins, const string &name,
                    const bool sumw2, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name, sumw2),f));
        }


        SigmaK0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            auto label_im_ng = [] (const unsigned ngamma)
            {
                const string label = std_ext::formatter() << ngamma << "#gamma IM [MeV]";
                return label;
            };


            AddTH1("KinFit Probability",      "probability",             "",       probbins,   "KinFitProb", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_prob(), f.TaggW());
            });
            AddTH1("TreeFit Probability",      "probability",             "",       probbins,   "TreeFitProb", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.SIG_prob(), f.TaggW());
            });

            AddTH1("6#gamma IM",label_im_ng(6), "", IMbins,"IM_6g", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM6g(), f.TaggW());
            });

            AddTH1("6#gamma IM fitted",label_im_ng(6), "", IMbins,"IM_6g_fit", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM6g(), f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", Bins(300,500,1900),"CBESUM",false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });

            AddTH2("Resonance Search 1","m(p #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300, 1050, 1700),Bins(300, 200, 900),"ppi0_2pi0", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for (auto i = 0u; i < pions.size() ; ++i)
                {
                    const auto N    = pions.at(i) + proton;
                    LorentzVec pipi({0,0,0},0);
                    for (auto j = 0u; j < pions.size() ; ++j)
                        if ( j != i )
                            pipi += pions.at(j);

                    h->Fill(N.M(),pipi.M(),f.TaggW());
                }
            });

            AddTH2("Resonance Search 1a","m^{2}(p #pi^{0}) [MeV]","m^{2}(2 #pi^{0}) [MeV]",Bins(300, std_ext::sqr(1050), std_ext::sqr(1700)),
                                                                                           Bins(300, std_ext::sqr( 200), std_ext::sqr( 800)),
                   "ppi0_2pi0_sqr", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for (auto i = 0u; i < pions.size() ; ++i)
                {
                    const auto N    = pions.at(i) + proton;
                    LorentzVec pipi({0,0,0},0);
                    for (auto j = 0u; j < pions.size() ; ++j)
                        if ( j != i )
                            pipi += pions.at(j);

                    h->Fill(N.M2(),pipi.M2(),f.TaggW());
                }
            });

            AddTH1("Resonance Search 2","m(p #pi^{0}) [MeV]","",Bins(300,  900, 1900),"ppi0", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for(auto i = 0u ; i < 3 ; ++i)
                {
                    const auto N    = pions.at(i) + f.Tree.EMB_proton();
                    h->Fill(N.M(),f.TaggW());
                }
            });

            AddTH2("Resonance Search 3","m(2 #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300,  0, 1000),Bins(300,    0, 1000),"2pi0_2pi0", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const vector<pair<size_t,size_t>> combinations = { { 0 , 1 } , { 0 , 2 } , { 1 , 2 } };
                const auto pions = f.Tree.SIG_pions();

                for ( size_t i = 0 ; i < 3 ; ++i)
                    for ( size_t j = 0 ; j < 3 ; ++j)
                    {
                        if ( i == j )
                            continue;
                        const auto ppM2  =(pions.at(combinations.at(i).first) + pions.at(combinations.at(i).second)).M();
                        const auto ppM1  =(pions.at(combinations.at(j).first) + pions.at(combinations.at(j).second)).M();
                        h->Fill(ppM2,ppM1,f.TaggW());
                    }
            });

            AddTH1("Resonance Search 4","m(K^{0}_{S} - candidate) [MeV]","",Bins(1000, 0 , 1000),"p2pi0", false,
                   [] (TH1* h, const Fill_t& f)
            {

                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                const auto Msigma2 = std_ext::sqr(ParticleTypeDatabase::SigmaPlus.Mass());
                double bestM2Diff = std_ext::inf;
                LorentzVec k0Cand({0,0,0},0);
                for( auto i = 0u ; i < 3 ; ++i)
                {
                    const auto N  = pions.at(i) + f.Tree.EMB_proton();
                    const auto m2d = std_ext::abs_diff(N.M2(),Msigma2);
                    if ( m2d < bestM2Diff)
                    {
                        bestM2Diff = m2d;
                        k0Cand = pions.at((i+1) % 3) + pions.at( (i+2) % 3);
                    }
                }

                h->Fill(k0Cand.M(),f.TaggW());
            });

            AddTH1("MC-true for reconstructed events","Tagger channel","",taggerBins,"effrecon", true,
                   [] (TH1* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(), f.TaggW());
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

        struct TreeCuts {

            static bool TreeFitProb(const Fill_t& f, const double p) noexcept {
                return     f.Tree.EMB_prob >  p;
            }
            static bool proton_MM(const Fill_t& f) noexcept {
                const auto width = 180.0;
                const auto mmpm = f.Tree.proton_MM();
                return (938.3 - width < mmpm && mmpm < 938.3 + width);
            }
            static bool dEk(const Fill_t& f, const double dEk)
            {
                return f.Tree.DiscardedEk() < dEk;
            }
        };

        static cuttree::Cuts_t<Fill_t> GetCuts() {

            using cuttree::MultiCut_t;

            cuttree::Cuts_t<Fill_t> cuts;

            const cuttree::Cut_t<Fill_t> ignore({"ignore", [](const Fill_t&){ return true; }});

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "IM 6g >  600 MeV", [](const Fill_t& f)
                                      {
                                          return f.Tree.SIG_IM6g() > 600;
                                      }
                                  }
                              });

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "discarded E_k == 0",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::dEk(f,10.0);
                                    }
                                  },
                                  { "discarded E_k < 20.",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::dEk(f,20.0);
                                    }
                                  },
                                  { "discarded E_k < 50.",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::dEk(f,50.0);
                                    }
                                  }
                              });


            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "prob > 0.01",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::TreeFitProb(f,0.01);
                                    }
                                  },
                                  { "prob > 0.05",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::TreeFitProb(f,0.05);
                                    }
                                  },
                                  { "prob > 0.1",
                                    [](const Fill_t& f)
                                    {
                                        return TreeCuts::TreeFitProb(f,0.1);
                                    }
                                  }
                              });



            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"Pi0PIDVeto==0",     [](const Fill_t& f) { return f.Tree.PionPIDVetoE() == 0;   }},
                                  {"Pi0PIDVeto<0.2",    [](const Fill_t& f) { return f.Tree.PionPIDVetoE() <  0.2; }}
                              });

            return cuts;
        }

    };

    plot::cuttree::Tree_t<MCTrue_Splitter<SigmaK0Hist_t>> signal_hists;

    static const string data_name;
    static const double binScale;
    unsigned nchannels;


    TTree* t = nullptr;
    WrapTree tree;

    sigmaPlus::SeenTree seenTree;
    TTree* seen;
    TH1D* hist_seenMC =  nullptr;


    sigmaPlus::RecTree recTree;
    TTree* rec;

    TH1D* eff;

    static const BinSettings taggerBins;


    virtual long long GetNumEntries() const override {return t->GetEntries();}

public:

    sigmaPlus_Plot(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        Plotter(name,input,opts),
        nchannels(ExpConfig::Setup::GetDetector<TaggerDetector_t>()->GetNChannels())
    {
        global_opts = opts;

        // load Main tree
        if(!input.GetObject(WrapTree::treeAccessName(),t))
            throw Exception("Input TTree not found");
        if(!tree.Matches(t))
            throw std::runtime_error("Tree branches don't match");
        tree.LinkBranches(t);


        // load efficiency related trees
        if(!input.GetObject(sigmaPlus::SeenTree::treeAccessName(),seen))
            throw Exception("Input TTree for seen mc not found");
        if(!seenTree.Matches(seen))
            throw std::runtime_error("Tree branches don't match");
        seenTree.LinkBranches(seen);

        if(!input.GetObject(sigmaPlus::RecTree::treeAccessName(),rec))
            throw Exception("Input TTree for reconstructed mc not found");
        if(!recTree.Matches(rec))
            throw std::runtime_error("Tree branches don't match");
        recTree.LinkBranches(rec);
        if (rec->GetEntries() != t->GetEntries())
            throw std::runtime_error("Tree size of main tree and rec tree dont match.");

        hist_seenMC = HistFac.makeTH1D("seenMC","Tagger channel","",taggerBins,"seenMC",true);
        for (long long en = 0 ; en < seen->GetEntries() ; ++en)
        {
            seen->GetEntry(en);
            hist_seenMC->Fill(seenTree.TaggerBin());
        }

        signal_hists = cuttree::Make<MCTrue_Splitter<SigmaK0Hist_t>>(HistFac);
    }


    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);
        rec->GetEntry(entry);

        cuttree::Fill<MCTrue_Splitter<SigmaK0Hist_t>>(signal_hists, {tree, recTree});
    }

    virtual void Finish() override{}
    virtual void ShowResult() override{}

    virtual ~sigmaPlus_Plot(){}

};




const BinSettings sigmaPlus_Plot::taggerBins(47); // TODO: Get from Setup
const double sigmaPlus_Plot::binScale  = 1.0;
const string sigmaPlus_Plot::data_name = "Data";

class sigmaPlus_FinalPlot: public Plotter {

protected:

    template<typename Hist_t>
    struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

        // Hist_t should have that type defined
        using Fill_t = typename Hist_t::Fill_t;



        MCTrue_Splitter(const HistogramFactory& histFac,
                        const cuttree::TreeInfo_t& treeInfo) :
            cuttree::StackedHists_t<Hist_t>(histFac, treeInfo)

        {
            using histstyle::Mod_t;


            this->GetHist(sigmaPlus::settings_t::Index_Data,
                          data_name, Mod_t::MakeDataPoints(kBlack));
            this->GetHist(sigmaPlus::settings_t::Index_Signal,
                          "Sig",  Mod_t::MakeLine(kRed, 2));
            this->GetHist(sigmaPlus::settings_t::Index_MainBkg,
                          "MainBkg",  Mod_t::MakeLine(kGreen, 2));


            this->GetHist(sigmaPlus::settings_t::Index_SumMC,
                          "Sum_MC", Mod_t::MakeLine(kBlack, 1));
            this->GetHist(sigmaPlus::settings_t::Index_BkgMC,
                          "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

            this->GetHist(sigmaPlus::settings_t::Index_brokenTree,
                          "brokenTree", Mod_t::MakeLine(kGray,1,-1));
            this->GetHist(sigmaPlus::settings_t::Index_unregTree,
                          "untaggedTree", Mod_t::MakeLine(kGray,1,-1));



            for (auto i_offset = 0u ; i_offset < sigmaPlus::otherBackgrounds.size() ; ++i_offset)
            {
                const auto mctrue = sigmaPlus::settings_t::Index_Offset + i_offset;
                this->GetHist(mctrue,
                              physics::sigmaPlus::getOtherChannelNames(mctrue),
                              Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-10), 1, kGray+1) );
            }

        }

        void Fill(const Fill_t& f) {

            const unsigned mctrue = f.Tree.MCTrue;
            const Hist_t& hist = this->GetHist(mctrue);
            hist.Fill(f);

            // handle MC_all and MC_bkg
            if(mctrue>0) {
                this->GetHist(3).Fill(f);
                if(mctrue >= 8 || mctrue == 2)
                    this->GetHist(4).Fill(f);
            }
        }
    };

    struct SigmaK0Hist_t {

        using Tree_t  = physics::sigmaPlus::PionProdTree;
        using RTree_t = sigmaPlus::RecTree;


        struct Fill_t {
            const Tree_t&  Tree;
            const RTree_t& RTree;
            Fill_t(const Tree_t& t, const RTree_t& rtree) : Tree(t), RTree(rtree) {}

            double TaggW() const {
                return Tree.Tagg_W;
            }

            vector<TLorentzVector> get2G(const vector<TSimpleParticle>& photons) const
            {
                vector<TLorentzVector> acc;
                const auto& permutation(Tree.SIG_combination());
                for (size_t i = 0 ; i < permutation.size(); i+=2)
                {
                    TLorentzVector gg(photons.at(permutation.at(i)));
                    gg += photons.at(permutation.at(i+1));
                    acc.emplace_back(gg);
                }
                return acc;
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

        const BinSettings taggerBins = sigmaPlus_FinalPlot::taggerBins;


        const BinSettings Ebins    = Bins(1000, 0, 1000);

        const BinSettings Chi2Bins = Bins(250, 0,   25);
        const BinSettings probbins = Bins(250, 0,   1);

        const BinSettings IMbins       = Bins(1000,  200, 1100);
        const BinSettings IMProtonBins = Bins(1000,  600, 1200);
        const BinSettings IM2g         = Bins(1000,    0,  360);

        const BinSettings pThetaBins = Bins( 200,  0,   80);
        const BinSettings pEbins     = Bins( 350,  0, 1200);

        HistogramFactory HistFac;

        void AddTH1(const string &title, const string &xlabel, const string &ylabel,
                    const BinSettings &bins, const string &name,
                    const bool sumw2, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name, sumw2),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel,
                    const BinSettings &xbins, const BinSettings& ybins, const string &name,
                    const bool sumw2, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name, sumw2),f));
        }


        SigmaK0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            auto label_im_ng = [] (const unsigned ngamma)
            {
                const string label = std_ext::formatter() << ngamma << "#gamma IM [MeV]";
                return label;
            };


            AddTH1("KinFit Probability",      "probability",             "",       probbins,   "KinFitProb", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_prob(), f.TaggW());
            });
            AddTH1("TreeFit Probability",      "probability",             "",       probbins,   "TreeFitProb", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.SIG_prob(), f.TaggW());
            });

            AddTH1("6#gamma IM",label_im_ng(6), "", IMbins,"IM_6g", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM6g(), f.TaggW());
            });

            AddTH1("6#gamma IM fitted",label_im_ng(6), "", IMbins,"IM_6g_fit", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM6g(), f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", Bins(300,500,1900),"CBESUM",false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });

            AddTH2("Resonance Search 1","m(p #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300, 1050, 1700),Bins(300, 200, 900),"ppi0_2pi0", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for (auto i = 0u; i < pions.size() ; ++i)
                {
                    const auto N    = pions.at(i) + proton;
                    LorentzVec pipi({0,0,0},0);
                    for (auto j = 0u; j < pions.size() ; ++j)
                        if ( j != i )
                            pipi += pions.at(j);

                    h->Fill(N.M(),pipi.M(),f.TaggW());
                }
            });

            AddTH2("Resonance Search 1a","m^{2}(p #pi^{0}) [MeV]","m^{2}(2 #pi^{0}) [MeV]",Bins(300, std_ext::sqr(1050), std_ext::sqr(1700)),
                                                                                           Bins(300, std_ext::sqr( 200), std_ext::sqr( 800)),
                   "ppi0_2pi0_sqr", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for (auto i = 0u; i < pions.size() ; ++i)
                {
                    const auto N    = pions.at(i) + proton;
                    LorentzVec pipi({0,0,0},0);
                    for (auto j = 0u; j < pions.size() ; ++j)
                        if ( j != i )
                            pipi += pions.at(j);

                    h->Fill(N.M2(),pipi.M2(),f.TaggW());
                }
            });

            AddTH1("Resonance Search 2","m(p #pi^{0}) [MeV]","",Bins(300,  900, 1900),"ppi0", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                for(auto i = 0u ; i < 3 ; ++i)
                {
                    const auto N    = pions.at(i) + f.Tree.EMB_proton();
                    h->Fill(N.M(),f.TaggW());
                }
            });

            AddTH2("Resonance Search 3","m(2 #pi^{0}) [MeV]","m(2 #pi^{0}) [MeV]",Bins(300,  0, 1000),Bins(300,    0, 1000),"2pi0_2pi0", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                const vector<pair<size_t,size_t>> combinations = { { 0 , 1 } , { 0 , 2 } , { 1 , 2 } };
                const auto pions = f.Tree.SIG_pions();

                for ( size_t i = 0 ; i < 3 ; ++i)
                    for ( size_t j = 0 ; j < 3 ; ++j)
                    {
                        if ( i == j )
                            continue;
                        const auto ppM2  =(pions.at(combinations.at(i).first) + pions.at(combinations.at(i).second)).M();
                        const auto ppM1  =(pions.at(combinations.at(j).first) + pions.at(combinations.at(j).second)).M();
                        h->Fill(ppM2,ppM1,f.TaggW());
                    }
            });

            AddTH1("Resonance Search 4","m(K^{0}_{S} - candidate) [MeV]","",Bins(1000, 0 , 1000),"p2pi0", false,
                   [] (TH1* h, const Fill_t& f)
            {

                const auto pions = f.Tree.SIG_pions();
                const auto proton = f.Tree.SIG_proton();

                const auto Msigma2 = std_ext::sqr(ParticleTypeDatabase::SigmaPlus.Mass());
                double bestM2Diff = std_ext::inf;
                LorentzVec k0Cand({0,0,0},0);
                for( auto i = 0u ; i < 3 ; ++i)
                {
                    const auto N  = pions.at(i) + f.Tree.EMB_proton();
                    const auto m2d = std_ext::abs_diff(N.M2(),Msigma2);
                    if ( m2d < bestM2Diff)
                    {
                        bestM2Diff = m2d;
                        k0Cand = pions.at((i+1) % 3) + pions.at( (i+2) % 3);
                    }
                }

                h->Fill(k0Cand.M(),f.TaggW());
            });

            AddTH1("MC-true for reconstructed events","Tagger channel","",taggerBins,"effrecon", true,
                   [] (TH1* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(), f.TaggW());
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

        struct TreeCuts {

            static bool TreeFitProb(const Fill_t& f, const double p) noexcept {
                return     f.Tree.EMB_prob >  p;
            }
            static bool TaggERange(const Fill_t& f,const IntervalD& taggEnInterval) noexcept {
                return taggEnInterval.Contains(f.Tree.Tagg_E);
            }
            static bool dEk(const Fill_t& f, const double dEk)
            {
                return f.Tree.DiscardedEk() < dEk;
            }
            static bool finalCuts(const Fill_t& f)
            {
                return f.Tree.SIG_IM6g() > 600 &&
                       TreeCuts::dEk(f,20.0) &&
                       TreeCuts::TreeFitProb(f,0.05) &&
                       f.Tree.PionPIDVetoE() == 0;
            }
        };

        static cuttree::Cuts_t<Fill_t> GetCuts() {

            using cuttree::MultiCut_t;

            cuttree::Cuts_t<Fill_t> cuts;
            const auto n_GammaBins = 10u;
            auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
            const auto eMax = tagger->GetPhotonEnergy(0);
            const auto eMin = tagger->GetPhotonEnergy( 47u - 1u);  // get from setup!!!
            const auto ebinWidth = (eMax - eMin )/ n_GammaBins;

            using cuttree::Cut_t;
            MultiCut_t<Fill_t> eGammaBins;
            for (auto bin = 0u ; bin < n_GammaBins ; ++bin) {
                const IntervalD egInterval(eMin + (bin * ebinWidth),eMin + ((bin+1) * ebinWidth));
                eGammaBins.emplace_back(Cut_t<Fill_t>{std_ext::formatter() << bin, [egInterval](const Fill_t& f)
                                                      {
                                                          return TreeCuts::finalCuts(f) &&
                                                                 TreeCuts::TaggERange(f,egInterval);
                                                      }});
            }
            cuts.push_back(eGammaBins);

            
            return cuts;
        }

    };

    plot::cuttree::Tree_t<MCTrue_Splitter<SigmaK0Hist_t>> signal_hists;

    static const string data_name;
    static const double binScale;


    TTree* t = nullptr;
    WrapTree tree;

    sigmaPlus::SeenTree seenTree;
    TTree* seen;
    TH1D* hist_seenMC =  nullptr;


    sigmaPlus::RecTree recTree;
    TTree* rec;

    TH1D* eff;

    static const BinSettings taggerBins;


    virtual long long GetNumEntries() const override {return t->GetEntries();}

public:

    sigmaPlus_FinalPlot(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        Plotter(name,input,opts)
    {
        global_opts = opts;

        // load Main tree
        if(!input.GetObject(WrapTree::treeAccessName(),t))
            throw Exception("Input TTree not found");
        if(!tree.Matches(t))
            throw std::runtime_error("Tree branches don't match");
        tree.LinkBranches(t);


        // load efficiency related trees
        if(!input.GetObject(sigmaPlus::SeenTree::treeAccessName(),seen))
            throw Exception("Input TTree for seen mc not found");
        if(!seenTree.Matches(seen))
            throw std::runtime_error("Tree branches don't match");
        seenTree.LinkBranches(seen);

        if(!input.GetObject(sigmaPlus::RecTree::treeAccessName(),rec))
            throw Exception("Input TTree for reconstructed mc not found");
        if(!recTree.Matches(rec))
            throw std::runtime_error("Tree branches don't match");
        recTree.LinkBranches(rec);
        if (rec->GetEntries() != t->GetEntries())
            throw std::runtime_error("Tree size of main tree and rec tree dont match.");

        hist_seenMC = HistFac.makeTH1D("seenMC","Tagger channel","",taggerBins,"seenMC",true);
        for (long long en = 0 ; en < seen->GetEntries() ; ++en)
        {
            seen->GetEntry(en);
            hist_seenMC->Fill(seenTree.TaggerBin());
        }

        signal_hists = cuttree::Make<MCTrue_Splitter<SigmaK0Hist_t>>(HistFac);
    }


    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);
        rec->GetEntry(entry);

        cuttree::Fill<MCTrue_Splitter<SigmaK0Hist_t>>(signal_hists, {tree, recTree});
    }

    virtual void Finish() override{}
    virtual void ShowResult() override{}

    virtual ~sigmaPlus_FinalPlot(){}
};

const BinSettings sigmaPlus_FinalPlot::taggerBins(47); // TODO: Get from Setup
const double sigmaPlus_FinalPlot::binScale  = 1.0;
const string sigmaPlus_FinalPlot::data_name = "Data";

AUTO_REGISTER_PLOTTER(sigmaPlus_Plot)
AUTO_REGISTER_PLOTTER(sigmaPlus_FinalPlot)
