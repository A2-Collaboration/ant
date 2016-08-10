#include "ProcessTaggEff.h"

#include "slowcontrol/SlowControlVariables.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

ProcessTaggEff::ProcessTaggEff(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    fillDebug(opts->Get<bool>("fillDebug", false))
{
    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if (!Tagger) throw std::runtime_error("No Tagger found");
    auto nchannels = Tagger->GetNChannels();

    slowcontrol::Variables::TaggerScalers->Request();
    byChannel = HistFac.makeTH1D("","channel","Frequency [Hz]",nchannels,"TaggerFreqs");
    mainTree.CreateBranches(HistFac.makeTTree("mainTree"));
    mainTree.TaggFreqs().resize(nchannels);

    if (fillDebug)
    {
        debugTree.CreateBranches(HistFac.makeTTree("debugTree"));
        debugTree.TaggFreqs().resize(nchannels);
    }

}

ProcessTaggEff::~ProcessTaggEff() {}

void ProcessTaggEff::ProcessEvent(const TEvent& ev, manager_t& )
{
    mainTree.LGFreq = 0; // LG not implemented yet.
    if ( seenEvents == 0 )
        mainTree.EvID = ev.Reconstructed().ID;

    if (fillDebug)
    {
        debugTree.LGFreq = 0;
        debugTree.EvID = ev.Reconstructed().ID;
    }

    const auto scalars = slowcontrol::Variables::TaggerScalers->Get();
    unsigned channel = 0;
    for (const auto& value: scalars)
    {
        byChannel->Fill(channel,value);

        mainTree.TaggFreqs().at(channel) += value;

        if (fillDebug) debugTree.TaggFreqs().at(channel) = value;

        channel++;
    }

    if (fillDebug)
        debugTree.Tree->Fill();

    seenEvents++;
}

void ProcessTaggEff::Finish()
{
    // normalize
    byChannel->Scale(1./seenEvents);
    for (auto& freq: mainTree.TaggFreqs())
        freq = freq / (1. * seenEvents);
    mainTree.LGFreq = mainTree.LGFreq / (1. * seenEvents);
    mainTree.Tree->Fill();


    LOG(INFO) << "Filled Tree for " << seenEvents << " events, starting at "<< std_ext::to_iso8601(mainTree.EvID().Timestamp) << "." << endl;
}

void ProcessTaggEff::ShowResult()
{
    byChannel->Draw();
}

AUTO_REGISTER_PHYSICS(ProcessTaggEff)
