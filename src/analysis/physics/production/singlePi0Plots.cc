#include "physics/scratch/wolfes/tools/tools.h"
#include "singlePi0.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"
#include "base/std_ext/string.h"

#include "base/Logger.h"

#include "TH1D.h"

#include "plot/CutTree.h"



using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

using namespace std;

using singlePi0_PlotBase = TreePlotterBase_t<singlePi0::PionProdTree>;

auto singlePi0Cut = [](const singlePi0::PionProdTree& tree)
{
    return (
                tree.Neutrals < 2
           );
};

// target density [1/mub]
const double targetDensity = 0.4249E6;


class singlePi0_Test: public singlePi0_PlotBase{

protected:
    TH1D* mPi0Before      = nullptr;
    TH1D* mPi0            = nullptr;

    TH1D* cutVar_Neutrals = nullptr;

    TH2D* countsraw       = nullptr;
    TH2D* countsCor       = nullptr;
    TH2D* xsec            = nullptr;

    TH2D* efficiencies    = nullptr;
    TH1D* taggerScalars   = nullptr;

    bool cut() const
    {
        return singlePi0Cut(tree);
    }

    unsigned nchannels;

    WrapTFileInput eff_input;



public:

    class AngularCollection_t {
    protected:
        const size_t NTaggerChannels;
        const BinSettings CosThetaBinning;
        const BinSettings IMBins;
        using histKey_t = pair<int,int>;
        using hists_t   = map<histKey_t,TH1D*>;
        hists_t _data;
    public:
        AngularCollection_t(
                const size_t nTaggerChannels,
                const BinSettings& thetaBinning,
                const BinSettings& imbins,
                const HistogramFactory& histFac):
            NTaggerChannels(nTaggerChannels),
            CosThetaBinning(thetaBinning),
            IMBins(imbins)
        {
            for (size_t ch = 0 ; ch < NTaggerChannels ; ++ch)
            {
                for ( size_t thetaBin = 0 ; thetaBin < CosThetaBinning.Bins() ; ++thetaBin)
                {
                    _data[{ch,thetaBin}] = histFac.makeTH1D(
                                                  "IM","im 2#gamma [MeV]","#", IMBins,
                                                  std_ext::formatter() << "im_" << ch << "_" << thetaBin);
                }
            }
        }

        int Fill(const size_t taggerChannel, const double cosTheta,
                 const double value, const double weight)
        {
            const auto thetaBin = CosThetaBinning.getBin(cosTheta);
            if (thetaBin == -1)
                return  -1;
            if (taggerChannel >= NTaggerChannels)
                return -1;

            return _data.at({taggerChannel,thetaBin})->Fill(value,weight);
        }

    };

    singlePi0_Test(const string& name, const WrapTFileInput& input,
                   OptionsPtr opts):
        singlePi0_PlotBase(name,input,opts),
        eff_input(opts->Get<string>("eff", ""))
    {
        auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
        if (!Tagger) throw std::runtime_error("No Tagger found");
        nchannels = Tagger->GetNChannels();



        LOG(INFO) << "Loading efficiencies for " << eff_input.FileNames() << ".";
        if(!eff_input.GetObject("singlePi0/eff2d",efficiencies))
            throw  std::runtime_error("Input TH1D for efficiencies not found");
        LOG(INFO) << "Loading scalar counts histogram";
        if(!input.GetObject("singlePi0/taggerScalars",taggerScalars))
            throw std::runtime_error("histogramm for taggerScalars not found");


        BinSettings taggerbins(nchannels);


        mPi0Before = HistFac.makeTH1D("before cuts","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));
        mPi0       = HistFac.makeTH1D("selection","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));


        cutVar_Neutrals = HistFac.makeTH1D("cut variable: # neutral candidates",
                                           "# neutrals","#",
                                           BinSettings(5));

        BinSettings costheta(efficiencies->GetNbinsY(),-1,1);
        const BinSettings egamma(efficiencies->GetNbinsX());

        countsraw = HistFac.makeTH2D("counts", "taggerChannel", "cos(#theta_{#pi^{0}})",
                                     taggerbins, costheta,
                                     "counts", true);
        countsCor = HistFac.makeTH2D("counts / ( #eta * l)",
                                     "taggerChannel","cos(#theta_{#pi^{0}})",
                                     taggerbins, costheta,
                                     "countsCor", true);
        xsec      = HistFac.makeTH2D("cross section", "taggerChannel", "cos(#theta_{#pi^{0}})",
                                     egamma, costheta,
                                     "xsec", true);
    }

    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);

        const auto taggW = tree.Tagg_W();

        mPi0Before->Fill(tree.IM2g(),taggW);

        cutVar_Neutrals->Fill(tree.Neutrals(),taggW);

        if (cut()) return;
        mPi0->Fill(tree.IM2g(),taggW);

        const auto ch = tree.Tagg_Ch();
        const auto effcorFac = tree.ExpLivetime() * tree.Tagg_Eff();
        const auto scalerCount = taggerScalars->GetBinContent(ch+1);

        if (effcorFac > 0 && scalerCount > 0)
        {
            countsraw->Fill(ch, tree.cosThetaPi0COMS(), taggW);
            countsCor->Fill(ch, tree.cosThetaPi0COMS(), taggW / effcorFac);
            xsec->Fill(ch,      tree.cosThetaPi0COMS(), taggW / ( scalerCount *effcorFac));
        }
    }

    virtual void Finish() override
    {
        xsec->Divide(efficiencies);
        xsec->Scale(targetDensity);
    }


    virtual void ShowResult() override
    {
        canvas("view")
                << mPi0Before
                << mPi0
                << cutVar_Neutrals
                << endc;
        canvas("cross sections")
                << drawoption("colz")
                << countsraw
                << countsCor
                << xsec
                << efficiencies
                << endc;
    }
};

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
//        TH3D* h3;
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

        void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, const bool sumw2, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name, sumw2),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, const bool sumw2, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name, sumw2),f));
        }
/*
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
*/

        SinglePi0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            AddTH1("TreeFit Probability", "probability", "", probbins, "TreeFitProb", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_prob(), f.TaggW());
            });

            AddTH1("2#gamma IM","2#gamma IM [MeV]", "", IM2g,"IM_2g", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM2g(), f.TaggW());
            });

            AddTH1("2#gamma IM fitted","2#gamma IM [MeV]", "", IM2g,"IM_2g_fit", true,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM2g(), f.TaggW());
            });

            AddTH1("MM proton","MM_{proton} [MeV]", "", IMProtonBins, "IM_p", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IMproton_MM(), f.TaggW());
            });

            AddTH1("DiscardedEk","E [MeV]", "#", DiscardedEkBins,"discEk", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.DiscardedEk(), f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", BinSettings(300,500,1900),"CBESUM", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH1("lifetime", "lifetime","", BinSettings(200,0,1),"lifetime", false,
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.ExpLivetime, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE", false,
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });


            AddTH1("#pi^0 - fitted", "cos(#theta)","#", cosThetaBins ,"costhetafit", false,
                   [] (TH1D* h, const Fill_t& f)
            {

                h->Fill(f.Tree.EMB_cosThetaPi0COMS(),f.TaggW());
            });

            AddTH2("reconstructed","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"recon_fit", true,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.Tagg_Ch(),f.Tree.EMB_cosThetaPi0COMS(),f.TaggW());
            });

            AddTH2("reconstructed - lifetime corrected","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"recon_cor", true,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.Tagg_Ch(),f.Tree.cosThetaPi0COMS(), f.TaggW() / 0.5744 ); // TODO fix lifetime = 0 issue!!!!!!!!!!!!!!! this is anaverage over all runfiles!
            });
            AddTH2("reconstructed - lifetime corrected - emb fitted","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"recon_cor_fit", true,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.Tagg_Ch(),f.Tree.EMB_cosThetaPi0COMS(), f.TaggW() / 0.5744 ); // TODO fix lifetime = 0 issue!!!!!!!!!!!!!!! this is anaverage over all runfiles!
            });

            AddTH2("eff_reconstructed_pi0","Tagger channel","cos(#theta_{#pi^{0}})",taggerBins, cosThetaBins,"effrecon_pi0", true,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(),f.RTree.CosThetaPi0(),f.TaggW());
            });

            AddTH2("eff_reconstructed","Tagger channel","cos(#theta_{lab})",taggerBins, cosThetaBins,"effreconcos", true,
                   []( TH2D* h, const Fill_t& f)
            {
                h->Fill(f.RTree.TaggerBin(), cos(f.RTree.Theta()),f.TaggW());
            });



//            AddTaggChVSthetaPlots();

        }

        void Fill(const Fill_t& f) const {
            h1.Fill(f);
            h2.Fill(f);

//            h3->Fill(
//                        f.Tree.EMB_cosThetaPi0COMS(),
//                        f.Tree.Tagg_Ch(),
//                        f.Tree.EMB_IM2g(),
//                        f.TaggW() / f.Tree.ExpLivetime()
//                    );
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
//            v.emplace_back(h3);
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
                                  { "dicardedEk<20",  [](const Fill_t& f) { return TreeCuts::DircardedEk(f, 20.);  }},
                                  { "dicardedEk<50",  [](const Fill_t& f) { return TreeCuts::DircardedEk(f, 50.);  }}
                              });

            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  { "EMB_prob>0.05", [](const Fill_t& f){ return TreeCuts::KinFitProb(f, 0.05); }},
                                  { "EMB_prob>0.1",  [](const Fill_t& f){ return TreeCuts::KinFitProb(f, 0.1);  }}
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"AllPhotonsInCB", [](const Fill_t& f) { return !TreeCuts::touchesHole(f); }},
                                  ignore
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"NoTouchesHole", [](const Fill_t& f) { return TreeCuts::allPhotonsInCB(f); }},
                                  ignore
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"Pi0PIDVeto==0",     [](const Fill_t& f) { return f.Tree.PionPIDVetoE() == 0;   }},
                                  {"Pi0PIDVeto<0.2",    [](const Fill_t& f) { return f.Tree.PionPIDVetoE() <  0.2; }}
                              });

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


        hist_seenMCcosTheta = HistFac.makeTH2D("seenMCcosTheta","Tagger channel","cos(#theta(#pi^{0}))", taggerBins, eff_cosThetaBins, "seenMCcosTheta");
        hist_seenMCTheta    = HistFac.makeTH2D("seenMCTheta",   "Tagger channel","#theta_{lab} [#circ]", taggerBins, eff_ThetaBins,    "seenMCTheta");
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
AUTO_REGISTER_PLOTTER(singlePi0_Test)
