#include "physics/scratch/wolfes/tools/tools.h"
#include "singlePi0.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"
#include "base/std_ext/string.h"

#include "base/Logger.h"

#include "TH1D.h"

#include "plot/CutTree.h"

#include "analysis/physics/scratch/wolfes/tools/tools.h"



using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

using namespace std;

using singlePi0_PlotBase = TreePlotterBase_t<singlePi0::PionProdTree>;


OptionsPtr global_opts = nullptr;
auto get_is_final = [](const OptionsPtr& opts) {return opts->Get<bool>("final", false);};
auto get_cuts_str = [](const OptionsPtr& opts) {return opts->Get<string>("cutstring","dicardedEk<20/EMB_prob>0.05/ignore/Pi0PIDVeto==0");};

using namespace ant::analysis::plot;
using WrapTree = singlePi0::PionProdTree;

class singlePi0_Plot: public Plotter{

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

            // TODO: derive this from channel map
            this->GetHist(0, "data", Mod_t::MakeDataPoints(kBlack));
            this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2));
            this->GetHist(2, "MainBkg",  Mod_t::MakeLine(kGreen, 2));
            // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
            // see also Fill()
            this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 1));
            this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+1, -1));

        }

        void Fill(const Fill_t& f) {

            const unsigned mctrue = f.Tree.MCTrue();

            using histstyle::Mod_t;

            auto get_bkg_name = [] (const unsigned) {
                return "unknown"; //(int(mctrue));
            };

            using histstyle::Mod_t;

            const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                            this->GetHist(mctrue,
                                                           get_bkg_name(mctrue),
                                                           Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-10), 1, kGray+1)
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

    struct SinglePi0Hist_t {

        using Tree_t  = singlePi0::PionProdTree;
        using RTree_t = singlePi0::RecTree;

        struct Fill_t {
            const Tree_t&  Tree;
            const RTree_t& RTree;

            Fill_t(const Tree_t& t, const RTree_t& rtree) :
                Tree(t), RTree(rtree){}

            double TaggW() const {
                return Tree.Tagg_W;
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

//        struct TaggChMgr : std::vector<TH2D*> {
//            using vector<TH2D*>::vector;
//            void Fill(const Fill_t& data) const
//            {
//                const auto tagg_channel = data.Tree.Tagg_Ch();
//                if (!singlePi0_Plot::taggerBins.Contains(tagg_channel))
//                    throw runtime_error(std_ext::formatter() << "Bad tagger channel: " << tagg_channel << " not in " << singlePi0_Plot::taggerBins );
//                this->at(tagg_channel)->Fill(data.Tree.EMB_cosThetaPi0COMS,data.Tree.EMB_IM2g,data.TaggW() / data.Tree.ExpLivetime);
//            }
//        };

        HistMgr<TH1D> h1;
        HistMgr<TH2D> h2;
        TH3D* h3 = nullptr;
//        TaggChMgr     taggChHists;

        const BinSettings probbins = BinSettings(250, 0,   1);
        const BinSettings taggerBins = singlePi0_Plot::taggerBins;

        const BinSettings IMbins       = BinSettings(1000,  200, 1100);
        const BinSettings IMProtonBins = BinSettings(1000,  600, 1200);
        const BinSettings IM2g         = BinSettings(1000,    0,  360);

        const BinSettings DiscardedEkBins = BinSettings(100);

        const BinSettings pThetaBins = BinSettings( 200,  0,   80);
        const BinSettings pEbins     = BinSettings( 350,  0, 1200);

        const BinSettings cosThetaBins   = singlePi0_Plot::eff_cosThetaBins;
        const BinSettings ThetaBins      = singlePi0_Plot::eff_ThetaBins;

        HistogramFactory HistFac;

        enum class addTo
        {
            final, overview, both
        };

        void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name,
                    const bool sumw2, const addTo add_to, fillfunc_t<TH1D> f) {
            auto is_final = get_is_final(global_opts);
            if ( (add_to == addTo::both) ||
                 (add_to == addTo::final && is_final) ||
                 (add_to == addTo::overview && !is_final))
            {
                h1.emplace_back(HistFiller_t<TH1D>(
                                    HistFac.makeTH1D(title, xlabel, ylabel, bins, name, sumw2),f));
            }
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name,
                    const bool sumw2, const addTo add_to, fillfunc_t<TH2D> f) {
            auto is_final = get_is_final(global_opts);
            if ( (add_to == addTo::both) ||
                 (add_to == addTo::final && is_final) ||
                 (add_to == addTo::overview && !is_final))
            {
                h2.emplace_back(HistFiller_t<TH2D>(
                                    HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name, sumw2),f));
            }
        }

        void AddTaggChVSthetaPlots()
        {
            h3 = HistFac.makeTH3D("im(2#gamma) life time corrected",
                                  "cos(#theta_{coms})",             "tagger channel",           "im(2#gamma) [MeV]",
                                  singlePi0_Plot::eff_cosThetaBins, singlePi0_Plot::taggerBins, IM2g,
                                  "finalPlot");
//            for ( auto i = 0u; i < singlePi0_Plot::taggerBins.Bins(); ++i)
//            {
//                taggChHists.emplace_back(HistFac.makeTH2D("IM 2#gamma lifetime corrected","cos(#theta_{coms})","IM 2#gamma [MeV]",
//                                                          singlePi0_Plot::eff_cosThetaBins,IM2g,
//                                                          std_ext::formatter() << "ch" << i));
//            }
        }


        SinglePi0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            AddTH1("TreeFit Probability", "probability", "", probbins, "TreeFitProb", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_prob(), f.TaggW());
            });

            AddTH1("2#gamma IM","2#gamma IM [MeV]", "", IM2g,"IM_2g", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM2g(), f.TaggW());
            });

            AddTH1("2#gamma IM fitted","2#gamma IM [MeV]", "", IM2g,"IM_2g_fit", true, addTo::both,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM2g(), f.TaggW());
            });

            AddTH1("MM proton","MM_{proton} [MeV]", "", IMProtonBins, "IM_p", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IMproton_MM(), f.TaggW());
            });

            AddTH1("DiscardedEk","E [MeV]", "#", DiscardedEkBins,"discEk", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.DiscardedEk(), f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", BinSettings(300,500,1900),"CBESUM", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH1("lifetime", "lifetime","", BinSettings(200,0,1),"lifetime", false, addTo::overview,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.ExpLivetime, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE", false, addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });


            AddTH1("#pi^0 - fitted", "cos(#theta)","#", cosThetaBins ,"costhetafit", true, addTo::both,
                   [] (TH1D* h, const Fill_t& f)
            {

                h->Fill(f.Tree.EMB_cosThetaPi0COMS(),f.TaggW());
            });

            AddTH1("#pi^0 - raw", "cos(#theta)","#", cosThetaBins ,"costheta", true, addTo::both,
                   [] (TH1D* h, const Fill_t& f)
            {

                h->Fill(f.Tree.cosThetaPi0COMS(),f.TaggW());
            });


            AddTH2("splitoffs","cos(#theta_{#pi^{0}})","# clusters",cosThetaBins, BinSettings(10),"splitoffs",false,addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.cosThetaPi0COMS(),f.Tree.NCands(),f.TaggW());
            } );

            AddTH2("reconstructed","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"recon", true, addTo::both,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.Tagg_Ch(),f.Tree.cosThetaPi0COMS(),f.TaggW());
            });

            AddTH2("reconstructed kin fitted","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"recon_fit", true, addTo::both,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.Tagg_Ch(),f.Tree.EMB_cosThetaPi0COMS(),f.TaggW());
            });


            AddTH2("eff_reconstructed_pi0","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"effrecon_pi0", true, addTo::both,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(),f.RTree.CosThetaPi0(),f.TaggW());
            });

            AddTH2("eff_reconstructed","Tagger channel","cos(#theta_{lab})",taggerBins, cosThetaBins,"effreconcos", true, addTo::both,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(), cos(f.RTree.Theta()),f.TaggW());
            });



            // ============================================    pulls   ===========================================================================================
            const BinSettings pullSettings(120,-10,10);
            AddTH2("pulls_photons","#theta_{lab} [#circ]","pulls",BinSettings(90,0,180),pullSettings,"pulls_photons",false, addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                for ( auto i = 0u ; i < f.Tree.photons().size() ; ++i)
                {
                    h->Fill(f.Tree.EMB_photons().at(i).Theta() * 180 / 3.14159 , f.Tree.EMB_pull_g_thetas().at(i),f.TaggW());
                }
            });
            AddTH2("pulls_photons_picoms","cos(#theta^{cms}_{#pi^{0}})","pulls",cosThetaBins,pullSettings,"pulls_photons_picoms",false, addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                for ( auto i = 0u ; i < f.Tree.photons().size() ; ++i)
                {
                    h->Fill(cos(f.Tree.cosThetaPi0COMS()), f.Tree.EMB_pull_g_thetas().at(i),f.TaggW());
                }
            });

            AddTH2("pulls_protons","#theta_{lab} [#circ]","pulls",BinSettings(90,0,180),pullSettings,"pulls_protons",false, addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().Theta() * 180 / 3.14159 , f.Tree.EMB_pull_p_theta(),f.TaggW());
            });
            AddTH2("pulls_protons_picoms","cos(#theta^{cms}_{#pi^{0}})","pulls",cosThetaBins,pullSettings,"pulls_protons_picoms",false, addTo::overview,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(cos(f.Tree.cosThetaPi0COMS()), f.Tree.EMB_pull_p_theta(),f.TaggW());
            });


            // ============================================    migration  ===========================================================================================
            AddTH2("pion migration lab sytem","mc #theta_{lab}","kin fit #theta_{lab}",BinSettings(180), BinSettings(180),"pionmig", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(std_ext::radian_to_degree(f.RTree.Theta()), std_ext::radian_to_degree(f.Tree.EMB_photonSum().Theta()),f.TaggW());
            });

            AddTH2("pion migration change","fit (#theta_{#pi^{0}})","#Delta(#theta_{#pi^{0}}) / (mc  #theta_{#pi^{0}})) ",BinSettings(180), BinSettings(100,-0.5,0.5),"pionmigdllab", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill( std_ext::radian_to_degree(f.Tree.EMB_photonSum().Theta()), (f.RTree.Theta() - f.Tree.EMB_photonSum().Theta()) / f.RTree.Theta() ,f.TaggW());
            });

            AddTH2("pion migration coms","mc cos(#theta_{#pi^{0}})","kin fit cos(#theta_{#pi^{0}})",cosThetaBins, cosThetaBins,"pionmigc", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.RTree.CosThetaPi0(), f.Tree.EMB_cosThetaPi0COMS() ,f.TaggW());
            });

            AddTH2("pion migration change","fit cos(#theta_{#pi^{0}})","#Deltacos(#theta_{#pi^{0}}) ",cosThetaBins, BinSettings(100,-0.1,0.1),"pionmigd", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill( f.Tree.EMB_cosThetaPi0COMS(), (f.RTree.CosThetaPi0() - f.Tree.EMB_cosThetaPi0COMS()) ,f.TaggW());
            });

            AddTH2("photons: #theta_{rec} - #theta_{mc}","cos(#theta_{#pi^{0}})","lab: fittet #theta - true #theta  {#gamma} [#circ]",BinSettings(100,-1,1),BinSettings(100,-25,25), "photonangles", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                const auto& fittedPhotons = f.Tree.EMB_photons();
                const auto& trueThetas   = f.RTree.gThetas();

                if ( fittedPhotons.size() != trueThetas.size())
                    return;
                if (trueThetas.size() != 2)
                    return;

                vector<double> fittedThetas(trueThetas.size());
                transform(fittedPhotons.begin(),fittedPhotons.end(), fittedThetas.begin(),
                          [](const TLorentzVector& g) { return g.Theta();});

                pair<bool,bool> smallestDiff;
                auto diff = std_ext::inf;

                assert(fittedThetas.size() == 2 && trueThetas.size() == 2);
                for (const auto ir: {false,true})
                    for (const auto it: {false,true})
                    {
                        const auto tdiff =( fittedThetas.at(ir) - trueThetas.at(it));
                        if (abs(tdiff) < abs(diff))
                        {
                            diff = tdiff;
                            smallestDiff = {ir,it};
                        }
                    }

                h->Fill(
                            f.RTree.CosThetaPi0(),
                            std_ext::radian_to_degree(diff),
                            f.TaggW()
                            );
                h->Fill(
                            f.RTree.CosThetaPi0(),
                            std_ext::radian_to_degree(fittedThetas.at(!(smallestDiff.first)) - trueThetas.at(!(smallestDiff.second))),
                            f.TaggW()
                            );
            });

            AddTH2("lbtocosthetagamma","lab #theta_{#gamma}","cos(#theta_{#pi^{0}})",BinSettings(100,0,180),BinSettings(100,-1,1), "convangletocosgamma", true, addTo::overview,
                   []( TH2D* h, const Fill_t& f)
            {
                for (const auto& g: f.RTree.gThetas())
                    h->Fill( std_ext::radian_to_degree(g), cos(f.RTree.Theta()),f.TaggW());
            });


            if (get_is_final(global_opts))
                AddTaggChVSthetaPlots();

        }

        void Fill(const Fill_t& f) const {
            h1.Fill(f);
            h2.Fill(f);

            if (get_is_final(global_opts))
            {
                h3->Fill(
                            f.Tree.EMB_cosThetaPi0COMS(),
                            f.Tree.Tagg_Ch(),
                            f.Tree.EMB_IM2g(),
                            f.TaggW()
                            );
            }
//            taggChHists.Fill(f);
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
            if (get_is_final(global_opts))
            {
                v.emplace_back(h3);
            }
//            for (auto& e: taggChHists)
//            {
//                v.emplace_back(e);
//            }
            return v;
        }


        struct TreeCuts {

            static bool KinFitProb(const Fill_t& f, const double p) noexcept
            {
                return     f.Tree.EMB_prob >  p;
            }
            static bool DircardedEk(const Fill_t& f, const double e) noexcept
            {
                return     f.Tree.DiscardedEk < e;
            }
            static bool allPhotonsInCB(const Fill_t& f) {
                bool isInside = true;
                const auto startCB = 23.0;
                for (const auto& g: f.Tree.photons())
                {
                    isInside = std_ext::radian_to_degree(g.Theta()) > startCB;
                }
                return isInside;
            }
            static bool touchesHole(const Fill_t& f)
            {
                for (const auto& g: f.Tree.photons())
                    if (g.TouchesHole) return  true;

                return f.Tree.proton().TouchesHole;
            }

            static bool onlyRealNeutral(const Fill_t& f) noexcept {
                return f.Tree.Neutrals == 2;
            }
            static bool pidVetoE(const Fill_t& f, const double eThresh) noexcept
            {
                return f.Tree.PionPIDVetoE() > eThresh;
            }
        };


        static cuttree::Cuts_t<Fill_t> GetCuts() {

            using cuttree::MultiCut_t;

            cuttree::Cuts_t<Fill_t> cuts;

            const cuttree::Cut_t<Fill_t> ignore({"ignore", [](const Fill_t&){ return true; }});

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "dicardedEk==0",  [](const Fill_t& f) { return TreeCuts::DircardedEk(f, 10.);  }},    // cluster threshold: 12 MeV...
                                  { "dicardedEk<20",  [](const Fill_t& f) { return TreeCuts::DircardedEk(f, 20.);  }},
                                  { "dicardedEk<50",  [](const Fill_t& f) { return TreeCuts::DircardedEk(f, 50.);  }}
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "EMB_prob>0.01", [](const Fill_t& f){ return TreeCuts::KinFitProb(f, 0.01); }},
                                  { "EMB_prob>0.05", [](const Fill_t& f){ return TreeCuts::KinFitProb(f, 0.05); }},
                                  { "EMB_prob>0.10", [](const Fill_t& f){ return TreeCuts::KinFitProb(f, 0.1);  }}
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"Pi0PIDVeto==0",     [](const Fill_t& f) { return f.Tree.PionPIDVetoE() == 0;   }},
                                  {"Pi0PIDVeto<0.2",    [](const Fill_t& f) { return f.Tree.PionPIDVetoE() <  0.2; }}
                              });

            if (get_is_final(global_opts))
            {
                return cuttree::Cuts_t<Fill_t>({
                                                   { {"final", [] (const Fill_t& f){
                                                          return
                                                          TreeCuts::DircardedEk(f,global_opts->Get<double>("dEk", 20.) ) &&
                                                          TreeCuts::KinFitProb(f, global_opts->Get<double>("prob", 0.05)) &&
                                                          f.Tree.PionPIDVetoE() <= global_opts->Get<double>("veto", 0.0);
                                                      }
                                                     } }
                                               });
            }

            return cuts;
        }

    };

    plot::cuttree::Tree_t<MCTrue_Splitter<SinglePi0Hist_t>> signal_hists;


    TTree* t = nullptr;
    WrapTree tree;
    unsigned nchannels;

    singlePi0::SeenTree seenTree;
    TTree* seen;
    TH2D* hist_seenMCcosTheta =  nullptr;
    TH2D* hist_seenMCTheta    =  nullptr;


    singlePi0::RecTree recTree;
    TTree* rec;

    TH2D* eff;

    static const BinSettings eff_cosThetaBins;
    static const BinSettings eff_ThetaBins;
    static const BinSettings taggerBins;


    virtual long long GetNumEntries() const override {return t->GetEntries();}


    // Plotter interface
public:
    singlePi0_Plot(const string& name, const WrapTFileInput& input, OptionsPtr opts):
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
        if(!input.GetObject(singlePi0::SeenTree::treeAccessName(),seen))
            throw Exception("Input TTree for seen mc not found");
        if(!seenTree.Matches(seen))
            throw std::runtime_error("Tree branches don't match");
        seenTree.LinkBranches(seen);

        if(!input.GetObject(singlePi0::RecTree::treeAccessName(),rec))
            throw Exception("Input TTree for reconstructed mc not found");
        if(!recTree.Matches(rec))
            throw std::runtime_error("Tree branches don't match");
        recTree.LinkBranches(rec);
        if (rec->GetEntries() != t->GetEntries())
            throw std::runtime_error("Tree size of main tree and rec tree dont match.");


        hist_seenMCcosTheta = HistFac.makeTH2D("seenMCcosTheta","Tagger channel","cos(#theta(#pi^{0}))", taggerBins, eff_cosThetaBins, "seenMCcosTheta",true);
        hist_seenMCTheta    = HistFac.makeTH2D("seenMCTheta",   "Tagger channel","#theta_{lab} [#circ]", taggerBins, eff_ThetaBins,    "seenMCTheta",true);
        for (long long en = 0 ; en < seen->GetEntries() ; ++en)
        {
            seen->GetEntry(en);
            hist_seenMCcosTheta->Fill(seenTree.TaggerBin(),seenTree.CosThetaPi0());
            hist_seenMCTheta->Fill(seenTree.TaggerBin(),std_ext::radian_to_degree(seenTree.Theta()));

        }

        signal_hists = cuttree::Make<MCTrue_Splitter<SinglePi0Hist_t>>(HistFac);
    }



    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);
        rec->GetEntry(entry);
        cuttree::Fill<MCTrue_Splitter<SinglePi0Hist_t>>(signal_hists, {tree, recTree});
    }

    virtual void Finish() override
    {

    }



    virtual void ShowResult() override{}

    virtual ~singlePi0_Plot(){}

};

const BinSettings singlePi0_Plot::eff_cosThetaBins(32,-1,1);
const BinSettings singlePi0_Plot::eff_ThetaBins(32,0,180);
const BinSettings singlePi0_Plot::taggerBins(47); // TODO: Get from Setup



AUTO_REGISTER_PLOTTER(singlePi0_Plot)
