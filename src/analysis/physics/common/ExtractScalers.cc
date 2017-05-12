#include "ExtractScalers.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;



ExtractScalers::ExtractScalers(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    MediateOver(opts->Get<seconds_t>("MediateOver",1)),
    nchannels(47),
    TDCcounts(47),
    Medians(47)
{


    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Clocks->Request();
    slowcontrol::Variables::Trigger->Request();
    slowcontrol::Variables::Beam->Request();

    scalers.CreateBranches(HistFac.makeTTree(treeName()));

    scalers.TDCRates().resize(nchannels);
    scalers.TDCRateErrors().resize(nchannels);
    scalers.TaggRates().resize(nchannels);
    scalers.TaggRateErrors().resize(nchannels);
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
    }

    if (scalers.AbsTime > MediateOver)
    {
        calcRates();
        resetAll();
    }


}

void ExtractScalers::calcRates()
{

}

size_t ExtractScalers::nChannels()
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");

    return Tagger->GetNChannels();
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
    scalers.TaggRates()  = slowcontrol::Variables::TaggerScalers->GetRates();

    scalers.AbsTime += slowcontrol::Variables::Clocks->GetBeampolmonClock();

    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        scalers.TDCRates().at(ch) = ( 1.0e6 * TDCcounts.at(ch)
                                          / slowcontrol::Variables::Clocks->GetExpClock() );

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
