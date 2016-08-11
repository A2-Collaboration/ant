#include "PhotonFlux.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

list<Variable::ProcessorPtr> PhotonFlux::GetNeededProcessors()
{
    mode = mode_t::PbGlass;
    return {Processors::PbGlass};
}

double PhotonFlux::Get() const
{
    if(mode == mode_t::PbGlass) {
        return Processors::PbGlass->PbGlassAcqu.Get();
    }
    LOG(WARNING) << "Bug: No processor for PbGlass!";
    return std::numeric_limits<double>::quiet_NaN();
}

