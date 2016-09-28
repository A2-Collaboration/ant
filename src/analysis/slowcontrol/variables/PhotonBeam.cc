#include "PhotonBeam.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> PhotonBeam::GetNeededProcessors() const
{
    return {Processors::Beampolmon,Processors::Beam};
}

double PhotonBeam::GetPbGlass() const
{
    return Processors::Beampolmon->PbGlass.Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();
}

double PhotonBeam::GetIonChamber() const
{
    return Processors::Beam->IonChamber.Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();
}
