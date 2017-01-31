#include "MCClusterECorr.h"

#include "base/Logger.h"

#include "plot/root_draw.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCClusterECorr::MCClusterECorr(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void MCClusterECorr::ProcessEvent(const TEvent& event, manager_t&)
{

}

void MCClusterECorr::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(MCClusterECorr)