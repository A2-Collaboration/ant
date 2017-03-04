#include "physics/common/MCTrueOverview.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCTrueOverview::MCTrueOverview(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
}


void MCTrueOverview::ProcessEvent(const TEvent& event, manager_t&)
{

}

void MCTrueOverview::ShowResult()
{

}


AUTO_REGISTER_PHYSICS(MCTrueOverview)