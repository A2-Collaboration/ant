#include "physics/common/MCGunCheck.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "tree/TParticle.h"

#include "plot/HistStyle.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1D.h"
#include "TTree.h"

#include "utils/Combinatorics.h"

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
    if(!event.MCTrue().ParticleTree)
        return;

    const auto particles = [](const TEvent& event)
    {
        const auto gPT = event.MCTrue().ParticleTree->Daughters();
        vector<TParticlePtr> ret(gPT.size());
        transform(gPT.begin(),gPT.end(),ret.begin(),
                  [](const TParticleTree_t& en){return en->Get();});
        return ret;
    }(event);

    for (const auto& gP: particles)
    {
        t.names().emplace_back(gP->Type().Name());
        t.thetas().push_back(gP->Theta());
        t.phis().push_back(gP->Phi());
    }


    if (particles.size() > 1)
    {
        auto combs = utils::makeCombination(particles,2);
        do {
            t.openings().push_back(TParticle::CalcAngle(combs.at(0),
                                                        combs.at(1)));
        } while(combs.next());
    }

    t.fillAndReset();

}


void MCGunCheck::Finish()
{

}

void MCGunCheck::ShowResult()
{
    auto tree = t.Tree;

    canvas("check")
            << TTree_drawable(tree,"cos(openings)")
            << TTree_drawable(tree,"openings * 180 / 3.1415")
            << drawoption("colz")
            << TTree_drawable(tree,"cos(thetas):phis")
            << TTree_drawable(tree,"thetas * 180 / 3.1415 :phis * 180 / 3.1415")
            << endc;
}




AUTO_REGISTER_PHYSICS(MCGunCheck)
