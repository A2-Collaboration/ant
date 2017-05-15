#include "physics/scratch/wolfes/tools/tools.h"
#include "singlePi0.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

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


class singlePi0_Efficiency: public DetectionEffciencyBase_t<singlePi0::PionProdTree>{

public:
    singlePi0_Efficiency(const string& name, const WrapTFileInput& input,
                         OptionsPtr opts):
        DetectionEffciencyBase_t<singlePi0::PionProdTree>(name,input,opts){}

    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);

        if (singlePi0Cut(tree)) return;

        efficiencies->Fill(tree.Tagg_Ch(),tree.Tagg_W());

    }

    virtual void Finish() override
    {
        efficiencies->Divide(seenMC);
    }

    virtual void ShowResult() override
    {
        canvas("efficiencies") << efficiencies << endc;
    }
};

class singlePi0_Test: public singlePi0_PlotBase{

protected:
    TH1D* mPi0Before      = nullptr;
    TH1D* mPi0            = nullptr;

    TH1D* cutVar_Neutrals = nullptr;

    TH1D* countsraw       = nullptr;
    TH1D* countsCor       = nullptr;
    TH1D* xsec            = nullptr;

    TH1D* efficiencies    = nullptr;
    TH1D* taggerScalars   = nullptr;

    bool cut() const
    {
        return singlePi0Cut(tree);
    }

    unsigned nchannels;

    WrapTFileInput eff_input;


    // Plotter interface
public:
    singlePi0_Test(const string& name, const WrapTFileInput& input,
                   OptionsPtr opts):
        singlePi0_PlotBase(name,input,opts),
        eff_input(opts->Get<string>("eff", ""))
    {
        auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
        if (!Tagger) throw std::runtime_error("No Tagger found");
        nchannels = Tagger->GetNChannels();



        LOG(INFO) << "Loading efficiencies for " << eff_input.FileNames() << ".";
        if(!eff_input.GetObject("singlePi0_Efficiency/eff",efficiencies))
            throw  std::runtime_error("Input TH1D for efficiencies not found");
        LOG(INFO) << "Loading scalar counts histogram";
        if(!input.GetObject("singlePi0/taggerScalars",taggerScalars))
            throw std::runtime_error("histogramm for taggerScalars not found");


        //        counts.resize(nchannels);

        mPi0Before = HistFac.makeTH1D("before cuts","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));
        mPi0       = HistFac.makeTH1D("selection","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));


        cutVar_Neutrals = HistFac.makeTH1D("cut variable: # neutral candidates",
                                           "# neutrals","#",
                                           BinSettings(5));

        countsraw = HistFac.makeTH1D("counts",
                                     "taggerChannel","# pi0 evts.",
                                     BinSettings(nchannels),
                                     "counts", true);
        countsCor = HistFac.makeTH1D("counts / ( #eta * l)",
                                     "taggerChannel","",
                                     BinSettings(nchannels),
                                     "countsCor", true);
        xsec      = HistFac.makeTH1D("cross section",
                                     "taggerChannel","cross section [mub]",
                                     BinSettings(nchannels),
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

        if (effcorFac > 0)
        {
            countsraw->Fill(ch, taggW);
            countsCor->Fill(ch, taggW / effcorFac);
            xsec->Fill(ch,      taggW / effcorFac);
        }
    }

    virtual void Finish() override
    {
        xsec->Divide(efficiencies);
        xsec->Divide(taggerScalars);
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
                << countsraw
                << countsCor
                << xsec
                << efficiencies
                << endc;
    }
};

using namespace ant::analysis::plot;

class singlePi0_Plot: public singlePi0_PlotBase{


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

            const unsigned mctrue = f.Tree.MCTrue;

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

        using Tree_t = singlePi0::PionProdTree;

        struct Fill_t {
            const Tree_t& Tree;

            Fill_t(const Tree_t& t) : Tree(t) {}

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

        HistMgr<TH1D> h1;
        HistMgr<TH2D> h2;

        const BinSettings probbins = BinSettings(250, 0,   1);

        const BinSettings IMbins       = BinSettings(1000,  200, 1100);
        const BinSettings IMProtonBins = BinSettings(1000,  600, 1200);
        const BinSettings IM2g         = BinSettings(1000,    0,  360);

        const BinSettings pThetaBins = BinSettings( 200,  0,   80);
        const BinSettings pEbins     = BinSettings( 350,  0, 1200);

        const BinSettings cosTheta   = BinSettings(30,-1,1);

        HistogramFactory HistFac;

        void AddTH1(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name, fillfunc_t<TH1D> f) {
            h1.emplace_back(HistFiller_t<TH1D>(
                                HistFac.makeTH1D(title, xlabel, ylabel, bins, name),f));
        }

        void AddTH2(const string &title, const string &xlabel, const string &ylabel, const BinSettings &xbins, const BinSettings& ybins, const string &name, fillfunc_t<TH2D> f) {
            h2.emplace_back(HistFiller_t<TH2D>(
                                HistFac.makeTH2D(title, xlabel, ylabel, xbins, ybins, name),f));
        }

        SinglePi0Hist_t(const HistogramFactory& hf, cuttree::TreeInfo_t): HistFac(hf)
        {
            AddTH1("TreeFit Probability",      "probability",             "",       probbins,   "TreeFitProb",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_prob, f.TaggW());
            });

            AddTH1("2#gamma IM","2#gamma IM [MeV]", "", IM2g,"IM_2g",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.IM2g, f.TaggW());
            });

            AddTH1("2#gamma IM fitted","2#gamma IM [MeV]", "", IM2g,"IM_2g_fit",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_IM2g, f.TaggW());
            });

            AddTH1("MM proton","MM_{proton} [MeV]", "", IMProtonBins, "IM_p",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.proton_MM().M(), f.TaggW());
            });



            AddTH1("Proton_MM_Angle", "Angle [#circ]","", BinSettings(200,0,40),"MM_pAngle",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.pMM_angle,f.TaggW());
            });

            AddTH1("CB_ESum", "EsumCB [MeV]","", BinSettings(300,500,1900),"CBESUM",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.CBESum, f.TaggW());
            });

            AddTH2("Fitted Proton","E^{kin}_{p} [MeV]","#theta_{p} [#circ]",pEbins,pThetaBins,"pThetaVsE",
                   [] (TH2D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.EMB_proton().E() - ParticleTypeDatabase::Proton.Mass(), std_ext::radian_to_degree(f.Tree.EMB_proton().Theta()), f.TaggW());
            });

            AddTH1("#pi^0", "cos(#theta)","#", cosTheta ,"costheta",
                   [] (TH1D* h, const Fill_t& f)
            {
                h->Fill(f.Tree.photonSum().CosTheta(),f.TaggW());
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

            static bool KinFitProb(const Fill_t& f) noexcept {
                return     f.Tree.EMB_prob >  0.1;
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
            static bool onlyRealNeutral(const Fill_t& f) noexcept {
                return f.Tree.Neutrals == 2;
            }
        };

        // Sig and Ref channel share some cuts...
        static cuttree::Cuts_t<Fill_t> GetCuts() {

            using cuttree::MultiCut_t;

            cuttree::Cuts_t<Fill_t> cuts;

            const cuttree::Cut_t<Fill_t> ignore({"ignore", [](const Fill_t&){ return true; }});


            cuts.emplace_back(MultiCut_t<Fill_t>{
                                 { "EMB_prob > 0.1", [](const Fill_t& f)
                                   {
                                       return TreeCuts::KinFitProb(f);
                                   }
                                 }
                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"all photons in CB", [](const Fill_t& f)
                                   {
                                       return TreeCuts::allPhotonsInCB(f);
                                   }
                                  },
                                  {"all #gamma neutral", [](const Fill_t& f)
                                   {
                                       return TreeCuts::onlyRealNeutral(f);
                                   }
                                  }

                              });
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"all photons in CB", [](const Fill_t& f)
                                   {
                                       return TreeCuts::allPhotonsInCB(f);
                                   }
                                  },
                                  {"all #gamma neutral", [](const Fill_t& f)
                                   {
                                       return TreeCuts::onlyRealNeutral(f);
                                   }
                                  }
                              });

             return cuts;
        }

    };

    plot::cuttree::Tree_t<MCTrue_Splitter<SinglePi0Hist_t>> signal_hists;

    // Plotter interface
public:

    singlePi0_Plot(const string& name, const WrapTFileInput& input, OptionsPtr opts):
        singlePi0_PlotBase(name,input,opts)
    {
        signal_hists = cuttree::Make<MCTrue_Splitter<SinglePi0Hist_t>>(HistFac);
    }


    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);
        cuttree::Fill<MCTrue_Splitter<SinglePi0Hist_t>>(signal_hists, {tree});
    }

    virtual void Finish() override{}
    virtual void ShowResult() override{}

    virtual ~singlePi0_Plot(){}

};

AUTO_REGISTER_PLOTTER(singlePi0_Efficiency)
AUTO_REGISTER_PLOTTER(singlePi0_Plot)
AUTO_REGISTER_PLOTTER(singlePi0_Test)
