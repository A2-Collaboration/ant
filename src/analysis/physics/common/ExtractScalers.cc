#include "ExtractScalers.h"

#include "slowcontrol/SlowControlVariables.h"
#include "plot/root_draw.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

ExtractScalers::ExtractScalers(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    MediateOver(opts->Get<seconds_t>("MediateOver",1)),
    histFac("hfac")
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");

    nchannels = Tagger->GetNChannels();

    auto bs = BinSettings(nchannels);
    hist_scalers = histFac.makeTH1D("scalars - e^{-} counts",    "channel no.","# per scaler block", bs);
    hist_tdchits = histFac.makeTH1D("tdc     - #gamma counts",      "channel no.","# per scaler block", bs);

    hist_scalers_rate = histFac.makeTH1D("scalars - e^{-} rate",    "channel no.","freq [Hz]", bs);
    hist_tdchits_rate = histFac.makeTH1D("tdc     - #gamma rate",    "channel no.","freq [Hz]", bs);


    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Clocks->Request();
    slowcontrol::Variables::Trigger->Request();
    slowcontrol::Variables::PhotonBeam->Request();

    scalers.CreateBranches(HistFac.makeTTree(treeName()));

    scalers.TDCRates().resize(nchannels);
    scalers.TaggRates().resize(nchannels);
    TDCcounts.resize(nchannels);
    scalers.AbsTime = 0;
}


ExtractScalers::~ExtractScalers() {}

void ExtractScalers::ProcessEvent(const TEvent& ev, manager_t& )
{
    scalers.nEvtsPerRead++;
    seenEvents++;

    processTaggerHits(ev);

    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        seenScalerBlocks++;

        processBlock();

        scalers.Tree->Fill();
        resetAll();
    }
}



void ExtractScalers::processTaggerHits(const TEvent &ev)
{
    for (const auto& taggerhit: ev.Reconstructed().TaggerHits)
    {
        TDCcounts.at(taggerhit.Channel)++;
    }
}

void ExtractScalers::processBlock()
{
    scalers.TaggRates()  = slowcontrol::Variables::TaggerScalers->Get();

    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        scalers.TDCRates().at(ch) = ( 1.0e6 * TDCcounts.at(ch)
                                          / slowcontrol::Variables::Clocks->GetExpClock() );
        hist_scalers_rate->Fill(ch,scalers.TaggRates().at(ch));
        hist_tdchits_rate->Fill(ch,scalers.TDCRates().at(ch));
    }

    clock = slowcontrol::Variables::Clocks->GetExpClock();
}

void ExtractScalers::resetAll()
{
    scalers.nEvtsPerRead = 0;
    for (auto& tdcC: TDCcounts)
        tdcC = 0;
}



void ExtractScalers::Finish()
{

}

void ExtractScalers::ShowResult()
{

}


AUTO_REGISTER_PHYSICS(ExtractScalers)
