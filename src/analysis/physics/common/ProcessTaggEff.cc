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
        slowcontrol::Variables::TaggerScalers->Request();
        byChannel = HistFac.makeTH1D("Tagger Scalars","channel","#",48);
        mainTree.CreateBranches(HistFac.makeTTree("mainTree"));
}

ProcessTaggEff::~ProcessTaggEff() {}

void ProcessTaggEff::ProcessEvent(const TEvent& ev, manager_t& )
{
    mainTree.AvgLGFreq = 0; // LG not implemented yet.
    if ( seenEvents == 0 ) mainTree.StartID = ev.Reconstructed().ID;

    const auto scalars = slowcontrol::Variables::TaggerScalers->Get();
    unsigned channel = 0;
    for (const auto& value: scalars)
    {
        byChannel->Fill(channel,value);
        mainTree.AvgTaggFreqs().at(channel) = value;
        channel++;
    }
    seenEvents++;
}

void ProcessTaggEff::Finish()
{
    // normalize
    byChannel->Scale(1./seenEvents);
    for (auto& freq: mainTree.AvgTaggFreqs())
        freq = freq / (1. * seenEvents);
    mainTree.AvgLGFreq = mainTree.AvgLGFreq / (1. * seenEvents);

    mainTree.Tree->Fill();

    LOG(INFO) << "Seen " << seenEvents << " Events";
}

void ProcessTaggEff::ShowResult()
{
    byChannel->Draw();
}

AUTO_REGISTER_PHYSICS(ProcessTaggEff)
