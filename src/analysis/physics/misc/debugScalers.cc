#include "debugScalers.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

debugScaler::debugScaler(const std::string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    auto nchannels = Tagger->GetNChannels();

    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Clocks->Request();
    slowcontrol::Variables::Trigger->Request();
    slowcontrol::Variables::Beam->Request();


    scalerReads.CreateBranches(HistFac.makeTTree("scalerReads"));

    scalerReads.TaggRates().resize(nchannels);
    scalerReads.TDCHits().resize(nchannels);
    scalerReads.CoincidentTDCHits().resize(nchannels);
    scalerReads.TaggTimings().resize(nchannels);
}


debugScaler::~debugScaler() {}

void debugScaler::ProcessEvent(const TEvent& ev, manager_t& )
{
    scalerReads.nEvtsPerRead++;
    seenEvents++;

    processTaggerHits(ev);

    if(slowcontrol::Variables::TaggerScalers->HasChanged())
    {
        processBlock(ev);
        scalerReads.Tree->Fill();
        LOG(INFO) << "ScalerRead ==>  "
                  << "n = " << scalerReads.nEvtsPerRead << ", "
                  << "N = " << scalerReads.LastID().Lower;
        scalerReads.nEvtsPerRead = 0;
    }
}

void debugScaler::Finish()
{
    LOG(INFO) << "Filled Tree for " << seenEvents << "." << endl; // << " events, starting at "<< std_ext::to_iso8601(mainTree.EvID().Timestamp) << "." << endl;
}

void debugScaler::ShowResult()
{
    canvas("check") << TTree_drawable(scalerReads.Tree, "PbRate")
                    << TTree_drawable(scalerReads.Tree, "ExpLivetime")
                    << TTree_drawable(scalerReads.Tree, "ExpTriggerRate")
                    << TTree_drawable(scalerReads.Tree, "ExpTriggerRate/L1TriggerRate")
                    << TTree_drawable(scalerReads.Tree, "ExpTriggerRate / PbRate")
                    << TTree_drawable(scalerReads.Tree, "ExpLivetime * PbRate / ExpTriggerRate")
                    << endc;

    canvas("check2") << TTree_drawable(scalerReads.Tree, "PbRate")
                    << TTree_drawable(scalerReads.Tree, "ExpLivetime")
                    << TTree_drawable(scalerReads.Tree, "ExpTriggerRate")
                    << TTree_drawable(scalerReads.Tree, "ExpLivetime / ExpTriggerRate>>deadtime(50,0,0)")
                    << TTree_drawable(scalerReads.Tree, "ExpTriggerRate / ( 1 - ExpLivetime)>>calcPbrate(50,0,0)")
                    << TTree_drawable(scalerReads.Tree, "PbRate / ( 1 + PbRate * (ExpLivetime / ExpTriggerRate))>>calcTrigger(50,0,0)")
                    << endc;
}

void debugScaler::processBlock(const TEvent& ev)
{
    const auto Scalers = slowcontrol::Variables::TaggerScalers->GetRates();
    unsigned channel = 0;
    for (const auto& value: Scalers)
    {
        scalerReads.TaggRates().at(channel) = value;
        channel++;
    }

    scalerReads.Exp1MHz = slowcontrol::Variables::Clocks->GetExpClock();
    scalerReads.BeamPolMon1MHz = slowcontrol::Variables::Clocks->GetBeampolmonClock();
    scalerReads.ExpLivetime = slowcontrol::Variables::Trigger->GetExpLivetime();
    scalerReads.PbRate = slowcontrol::Variables::Beam->GetPbGlass();
    scalerReads.LastID = ev.Reconstructed().ID;
    scalerReads.ExpTriggerRate = slowcontrol::Variables::Trigger->GetExpTrigger();
    scalerReads.L1TriggerRate = slowcontrol::Variables::Trigger->GetL1Trigger();
}

void debugScaler::processTaggerHits(const TEvent &ev)
{
    for (const auto& taggerhit: ev.Reconstructed().TaggerHits)
    {
        scalerReads.TaggTimings().at(taggerhit.Channel).emplace_back(taggerhit.Time);
    }
}

AUTO_REGISTER_PHYSICS(debugScaler)
