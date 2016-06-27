#include "DebugTaggerScalars.h"

#include "slowcontrol/SlowControlVariables.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

DebugTaggerScalars::DebugTaggerScalars(const std::string& name, OptionsPtr opts) :
    Physics(name, opts)
{
        slowcontrol::Variables::TaggerScalers->Request();
        byChannel = HistFac.makeTH1D("Tagger Scalars","channel","#",48);
}

DebugTaggerScalars::~DebugTaggerScalars() {}

void DebugTaggerScalars::ProcessEvent(const TEvent& , manager_t& )
{
    const auto scalars = slowcontrol::Variables::TaggerScalers->Get();
    unsigned channel = 0;
    for (const auto& value: scalars)
    {
        byChannel->Fill(channel,value);
        channel++;
    }
    seenEvents++;
}

void DebugTaggerScalars::Finish()
{
    byChannel->Scale(1./seenEvents);
    LOG(INFO) << "Seen " << seenEvents << " Events";
}

void DebugTaggerScalars::ShowResult()
{
    byChannel->Draw();
}

AUTO_REGISTER_PHYSICS(DebugTaggerScalars)
