#include "physics/scratch/wolfes/tools/tools.h"
#include "singlePi0.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

#include "base/Logger.h"

#include "TH1D.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

using namespace std;

using singlePi0_PlotBase = TreePlotterBase_t<singlePi0::PionProdTree>;


class singlePi0_Test: public singlePi0_PlotBase{

protected:
    TH1D* mPi0Before      = nullptr;
    TH1D* mPi0            = nullptr;

    TH1D* cutVar_SIG_prob = nullptr;
    TH1D* cutVar_Neutrals = nullptr;

    TH1D* countsraw       = nullptr;
    TH1D* countsCor       = nullptr;
    TH1D* xsec            = nullptr;


    bool cut() const
    {
        return (
                    tree.SIG_prob < 0.1 &&
                    tree.Neutrals != 2
               );
    }

    unsigned nchannels;

//    vector<long long> counts;


    // Plotter interface
public:
    singlePi0_Test(const string& name, const WrapTFileInput& input,
                   OptionsPtr opts):
        singlePi0_PlotBase(name,input,opts)
    {
        auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
        if (!Tagger) throw std::runtime_error("No Tagger found");
        nchannels = Tagger->GetNChannels();
//        counts.resize(nchannels);

        mPi0Before = HistFac.makeTH1D("before cuts","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));
        mPi0       = HistFac.makeTH1D("selection","m(#pi^{0}) [MeV]","#",
                                      BinSettings(200,50,220));

        cutVar_SIG_prob  = HistFac.makeTH1D("cut variable: probability",
                                            "prob","#",
                                            BinSettings(200,0,1));
        cutVar_Neutrals = HistFac.makeTH1D("cut variable: # neutral candidates",
                                           "# neutrals","#",
                                           BinSettings(5));

        countsraw = HistFac.makeTH1D("counts",
                                     "taggerChannel","# pi0 evts.",
                                     BinSettings(nchannels),
                                     "counts");
        countsCor = HistFac.makeTH1D("counts / scaler_rate * exp lifetime",
                                     "taggerChannel","",
                                     BinSettings(nchannels),
                                     "countsCor",
                                     true);
        xsec      = HistFac.makeTH1D("counts / l",
                                     "taggerChannel","",
                                     BinSettings(nchannels),
                                     "xsex",
                                     true);
    }

    virtual void ProcessEntry(const long long entry) override
    {
        t->GetEntry(entry);

        mPi0Before->Fill(tree.IM2g());

        cutVar_Neutrals->Fill(tree.Neutrals());
        cutVar_SIG_prob->Fill(tree.SIG_prob());

        if (cut()) return;
        mPi0->Fill(tree.IM2g());

        const auto ch = tree.Tagg_Ch();
        const auto scRateLT = tree.TaggRates().at(ch) * tree.ExpLivetime();
        const auto lumi = scRateLT * tree.Tagg_Eff();//per channel taggeff!!!

        if (scRateLT != 0)
        {
            countsraw->Fill(ch);
            countsCor->Fill(ch,1/scRateLT);
            xsec->Fill(ch,1/lumi);
        }
    }


    virtual void ShowResult() override
    {
        canvas("view")
                << mPi0Before
                << mPi0
                << cutVar_Neutrals
                << cutVar_SIG_prob
                << endc;
        canvas("cross sections")
                << countsraw
                << countsCor
                << xsec
                << endc;
    }
};

AUTO_REGISTER_PLOTTER(singlePi0_Test)
