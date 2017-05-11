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

auto singlePi0Cut = [](const singlePi0::PionProdTree& tree)
{
    return (
                tree.SIG_prob < 0.1 &&
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

    TH1D* cutVar_SIG_prob = nullptr;
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

        cutVar_SIG_prob  = HistFac.makeTH1D("cut variable: probability",
                                            "prob","#",
                                            BinSettings(200,0,1));
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
        cutVar_SIG_prob->Fill(tree.SIG_prob(),taggW);

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
                << cutVar_SIG_prob
                << endc;
        canvas("cross sections")
                << countsraw
                << countsCor
                << xsec
                << efficiencies
                << endc;
    }
};

AUTO_REGISTER_PLOTTER(singlePi0_Efficiency)
AUTO_REGISTER_PLOTTER(singlePi0_Test)
