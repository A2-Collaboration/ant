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
    auto nchannels = Tagger->GetNChannels();

    slowcontrol::Variables::TaggerScalers->Request();
    slowcontrol::Variables::FreeRates->Request();
    taggerChannels = HistFac.makeTH1D("","channel","Frequency [Hz]",nchannels,"TaggerFreqs");

    scalarReads.CreateBranches(HistFac.makeTTree("scalarReads"));

    scalarReads.TaggRates().resize(nchannels);
    scalarReads.TDCHits().resize(nchannels);
    scalarReads.CoincidentTDCHits().resize(nchannels);
    scalarReads.TaggTimings().resize(nchannels);
}


ProcessTaggEff::~ProcessTaggEff() {}

void ProcessTaggEff::ProcessEvent(const TEvent& ev, manager_t& )
{
    scalarReads.nEvtsPerRead++;

    if(slowcontrol::Variables::FreeRates->HasChanged())
    {

        const auto scalars = slowcontrol::Variables::TaggerScalers->Get();
        unsigned channel = 0;
        for (const auto& value: scalars)
        {
            taggerChannels->Fill(channel,value * scalarReads.nEvtsPerRead);
            scalarReads.TaggRates().at(channel) = value;
            channel++;
        }

        scalarReads.LGRate += slowcontrol::Variables::FreeRates->GetPbGlass();
        scalarReads.LastID = ev.Reconstructed().ID;
        scalarReads.Tree->Fill();
        LOG(INFO) << "ScalarRead ==>  "
                  << "n = " << scalarReads.nEvtsPerRead << ", "
                  << "N = " << scalarReads.LastID().Lower;
        scalarReads.nEvtsPerRead = 0;
    }

    seenEvents++;
}

void ProcessTaggEff::Finish()
{
    LOG(INFO) << "Filled Tree for " << seenEvents << "." << endl; // << " events, starting at "<< std_ext::to_iso8601(mainTree.EvID().Timestamp) << "." << endl;
}

void ProcessTaggEff::ShowResult()
{
    taggerChannels->Draw();
}

AUTO_REGISTER_PHYSICS(ProcessTaggEff)
