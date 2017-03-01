#include "physics/common/MCGunCheck.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"


#include "plot/HistStyle.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1D.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCGunCheck::MCGunCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
    t.CreateBranches(HistFac.makeTTree("tree"));
}


void MCGunCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    LOG(INFO) << event.MCTrue().Candidates.size() << "cands";
}


void MCGunCheck::Finish()
{

}

void MCGunCheck::ShowResult()
{

}




AUTO_REGISTER_PHYSICS(MCGunCheck)
