#include "ProcessTaggEff.h"

#include "slowcontrol/SlowControlVariables.h"
#include "plot/root_draw.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

ProcessTaggEff::ProcessTaggEff(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
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
    slowcontrol::Variables::FreeRates->Request();

    scalerReads.CreateBranches(HistFac.makeTTree(treeName()));

    scalerReads.TDCRates().resize(nchannels);
    scalerReads.TDCCounts().resize(nchannels);
    scalerReads.TaggCounts().resize(nchannels);
    scalerReads.TaggRates().resize(nchannels);

    scalerReads.TaggTimings().resize(nchannels);
}


ProcessTaggEff::~ProcessTaggEff() {}

void ProcessTaggEff::ProcessEvent(const TEvent& ev, manager_t& )
{
    scalerReads.nEvtsPerRead++;
    seenEvents++;

    processTaggerHits(ev);

    if(slowcontrol::Variables::FreeRates->HasChanged())
    {
        scalerReads.EvID = ev.Reconstructed().ID;
        seenScalerBlocks++;

        processBlock();
        scalerReads.Tree->Fill();

        resetAll();
    }
}



void ProcessTaggEff::processTaggerHits(const TEvent &ev)
{
    for (const auto& taggerhit: ev.Reconstructed().TaggerHits)
    {
        scalerReads.TDCCounts().at(taggerhit.Channel)++;
        scalerReads.TaggTimings().at(taggerhit.Channel).emplace_back(taggerhit.Time);

        hist_tdchits->Fill(taggerhit.Channel);
    }
}

void ProcessTaggEff::processBlock()
{
    scalerReads.TaggRates()  = slowcontrol::Variables::TaggerScalers->Get();

    for ( auto ch = 0u ; ch < nchannels ; ++ch)
    {
        scalerReads.TaggCounts().at(ch) = slowcontrol::Variables::TaggerScalers->GetCounts().at(ch);
        scalerReads.TDCRates().at(ch) = ( 1.0e6 * scalerReads.TDCCounts().at(ch)
                                          / slowcontrol::Variables::FreeRates->GetExpClock() );
        hist_scalers->Fill(ch,scalerReads.TaggCounts().at(ch));
        hist_scalers_rate->Fill(ch,scalerReads.TaggRates().at(ch));
        hist_tdchits_rate->Fill(ch,scalerReads.TDCRates().at(ch));
    }

    scalerReads.Clock = slowcontrol::Variables::FreeRates->GetExpClock();
    scalerReads.ExpLivetime = slowcontrol::Variables::FreeRates->GetExpLivetime();
    scalerReads.PbRate = slowcontrol::Variables::FreeRates->GetPbGlass();
    scalerReads.ExpTriggerRate = slowcontrol::Variables::FreeRates->GetExpTrigger();
}

void ProcessTaggEff::resetAll()
{
    scalerReads.nEvtsPerRead = 0;

    for ( auto ch = 0u ; ch < nchannels ; ++ch )
    {
        scalerReads.TDCCounts().at(ch) = 0;
        scalerReads.TaggTimings().at(ch).clear();
    }
}



void ProcessTaggEff::Finish()
{
    hist_scalers_rate->Scale(1.0 / seenScalerBlocks);
    hist_tdchits_rate->Scale(1.0 / seenScalerBlocks);

    LOG(INFO) << "Filled Tree for " << seenEvents
              << "in " << seenScalerBlocks
              << " scalerBlocks." << endl;
}

void ProcessTaggEff::ShowResult()
{

    canvas("check") << TTree_drawable(scalerReads.Tree, "PbRate")
                    << TTree_drawable(scalerReads.Tree, "TDCRates")
                    <<  endc;

    hist_scalers->SetLineColor(kRed);
    hist_scalers_rate->SetLineColor(kRed);
    canvas("channels") << padoption::Legend << hist_scalers      << samepad << hist_tdchits
                       << padoption::Legend << hist_scalers_rate << samepad << hist_tdchits_rate
                       << endc;
}


AUTO_REGISTER_PHYSICS(ProcessTaggEff)
