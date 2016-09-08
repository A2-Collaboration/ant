#include "etaprime_sergey.h"

#include "plot/root_draw.h"
#include "base/std_ext/misc.h"
#include "base/Logger.h"
#include "utils/particle_tools.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

EtapSergey::EtapSergey(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void EtapSergey::ProcessEvent(const TEvent& event, manager_t&)
{
    const TEventData& data = event.Reconstructed();

    auto result = fitter_sergey.Process(data);

}

void EtapSergey::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(EtapSergey)
