#include "Beam.h"

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"

using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;


list<Variable::ProcessorPtr> Beam::GetNeededProcessors() const
{
    return {Processors::Beampolmon,Processors::Beam};
}

double Beam::GetPbGlass() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::Beampolmon->PbGlass.Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();
}

double Beam::GetIonChamber() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::Beam->IonChamber.Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();
}

double Beam::GetFaradayCup() const
{
    if(!slowcontrol_provided)
        return 1.0;

    return Processors::Beampolmon->FaradayCup.Get();
}
