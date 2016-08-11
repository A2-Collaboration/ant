#include "PhotonFlux.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

void PhotonFlux::Init()
{
    mode = mode_t::PbGlass;
}

list<Variable::ProcessorPtr> PhotonFlux::GetNeededProcessors() const
{
    return {Processors::Beampolmon};
}

double PhotonFlux::Get() const
{
    if(mode == mode_t::PbGlass) {
        return Processors::Beampolmon->PbGlass.Get();
    }
    LOG(WARNING) << "Bug: No processor for PbGlass!";
    return std::numeric_limits<double>::quiet_NaN();
}

