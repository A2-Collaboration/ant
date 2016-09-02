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
    auto nchannels = Tagger->GetNChannels();

    auto bs = BinSettings(nchannels);
    hist_scalers = histFac.makeTH1D("scalars - e^{-} counts",    "channel no.","# per scaler block", bs);
    hist_tdchits = histFac.makeTH1D("tdc     - #gamma counts",      "channel no.","# per scaler block", bs);

    hist_scalers_rate = histFac.makeTH1D("scalars - e^{-} rate",    "channel no.","freq [Hz]", bs);
    hist_tdchits_rate = histFac.makeTH1D("tdc     - #gamma rate",    "channel no.","freq [Hz]", bs);


    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::FreeRates->Request();

    scalerReads.CreateBranches(HistFac.makeTTree(treeName()));

    scalerReads.TaggRates().resize(nchannels);
    scalerReads.TDCHits().resize(nchannels);
//    scalerReads.CoincidentTDCHits().resize(nchannels);
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
        processBlock(ev);
        scalerReads.Tree->Fill();
        LOG(INFO) << "ScalerRead ==>  "
                  << "n = " << scalerReads.nEvtsPerRead << ", "
                  << "N = " << scalerReads.LastID().Lower;
        scalerReads.nEvtsPerRead = 0;
        for (auto& nhits: scalerReads.TDCHits())
            nhits = 0;
        for (auto& taggerhits: scalerReads.TaggTimings())
            taggerhits.clear();
    }
}

void ProcessTaggEff::Finish()
{
    LOG(INFO) << "Filled Tree for " << seenEvents << "." << endl; // << " events, starting at "<< std_ext::to_iso8601(mainTree.EvID().Timestamp) << "." << endl;
}

void ProcessTaggEff::ShowResult()
{

    canvas("check") << TTree_drawable(scalerReads.Tree, "PbRate")
                    << TTree_drawable(scalerReads.Tree, "TDCHits")
                    <<  endc;
    hist_scalers->SetLineColor(kRed);
    canvas("channels") << padoption::Legend << hist_scalers      << samepad << hist_tdchits
                       << padoption::Legend << hist_scalers_rate << samepad << hist_tdchits_rate
                       << endc;
}

void ProcessTaggEff::processBlock(const TEvent& ev)
{
    scalerReads.TaggRates() = slowcontrol::Variables::TaggerScalers->GetCounts();

    for ( auto ch = 0u ; ch < scalerReads.TaggRates().size() ; ++ch)
        hist_scalers->Fill(ch,scalerReads.TaggRates().at(ch));

    scalerReads.ExpLivetime = slowcontrol::Variables::FreeRates->GetExpLivetime();
    scalerReads.PbRate = slowcontrol::Variables::FreeRates->GetPbGlass();
    scalerReads.LastID = ev.Reconstructed().ID;
    scalerReads.ExpTriggerRate = slowcontrol::Variables::FreeRates->GetExpTrigger();
}

void ProcessTaggEff::processTaggerHits(const TEvent &ev)
{
    for (const auto& taggerhit: ev.Reconstructed().TaggerHits)
    {
        scalerReads.TDCHits().at(taggerhit.Channel)++;
        hist_tdchits->Fill(taggerhit.Channel);
        scalerReads.TaggTimings().at(taggerhit.Channel).emplace_back(taggerhit.Time);
    }
}

AUTO_REGISTER_PHYSICS(ProcessTaggEff)
