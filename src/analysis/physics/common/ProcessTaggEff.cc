#include "ProcessTaggEff.h"

#include "slowcontrol/SlowControlVariables.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

ProcessTaggEff::ProcessTaggEff(const std::string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");

    nchannels = Tagger->GetNChannels();

    auto bs = BinSettings(nchannels);
    hist_scalers = HistFac.makeTH1D("scalars - e^{-} counts",    "channel no.","# per scaler block", bs);
    hist_tdchits = HistFac.makeTH1D("tdc     - #gamma counts",      "channel no.","# per scaler block", bs);

    hist_scalers_rate = HistFac.makeTH1D("scalars - e^{-} rate",    "channel no.","freq [Hz]", bs);
    hist_tdchits_rate = HistFac.makeTH1D("tdc     - #gamma rate",    "channel no.","freq [Hz]", bs);

    hist_tdc_times = HistFac.makeTH1D("tdc time", "tdc time","Counts", BinSettings(600,-150,150));
    hist_tdc_times_ch = HistFac.makeTH2D("tdc time v channel", "tdc time ","Channel", BinSettings(600,-150,150), bs);



    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::Clocks->Request();
    slowcontrol::Variables::Trigger->Request();
    slowcontrol::Variables::Beam->Request();

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

    if(slowcontrol::Variables::TaggerScalers->HasChanged())
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

        hist_tdc_times->Fill(taggerhit.Time);
        hist_tdc_times_ch->Fill(taggerhit.Time, taggerhit.Channel);
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
                                          / slowcontrol::Variables::Clocks->GetExpClock() );
        hist_scalers->Fill(ch,scalerReads.TaggCounts().at(ch));
        hist_scalers_rate->Fill(ch,scalerReads.TaggRates().at(ch));
        hist_tdchits_rate->Fill(ch,scalerReads.TDCRates().at(ch));
    }

    scalerReads.Clock = slowcontrol::Variables::Clocks->GetExpClock();
    scalerReads.ExpLivetime = slowcontrol::Variables::Trigger->GetExpLivetime();
    scalerReads.PbRate = slowcontrol::Variables::Beam->GetPbGlass();
    scalerReads.ExpTriggerRate = slowcontrol::Variables::Trigger->GetExpTrigger();
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
              << " in " << seenScalerBlocks
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
