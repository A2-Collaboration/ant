#include "physics/common/MCPhotonPairCheck.h"

#include "plot/root_draw.h"
#include "utils/particle_tools.h"
#include "tree/TParticle.h"

#include "plot/HistStyle.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1D.h"
#include "TTree.h"

#include "expconfig/ExpConfig.h"

#include "utils/combinatorics.h"



#include <cmath>
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCPhotonPairCheck::MCPhotonPairCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts),
    a2geo()
{
    t.CreateBranches(HistFac.makeTTree("tree"));
}


void MCPhotonPairCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& particleTree   = event.MCTrue().ParticleTree;

    if(!particleTree)
        return;


    // get photons;
    vector<TParticlePtr> photons;
    particleTree->Map_nodes([&photons](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::Photon)
            photons.emplace_back(pt->Get());
    });


    for (const auto& gP: photons)
    {
        t.p().emplace_back(*gP);
        const auto detHit = a2geo.DetectorFromAngles(gP->Theta(),
                                                     gP->Phi()    );
        if (detHit == Detector_t::Type_t::CB)
            t.hitsCB++;
        if (detHit == Detector_t::Type_t::TAPS)
            t.hitsTAPS++;
    }
    t.multiplicity = photons.size();


    if (photons.size() > 1)
    {
        auto combs = utils::makeCombination(photons,2);
        do {
            t.openings().push_back(TParticle::CalcAngle(combs.at(0),
                                                        combs.at(1)));
        } while(combs.next());
    }



    t.fillAndReset();

}


void MCPhotonPairCheck::Finish()
{

}

void MCPhotonPairCheck::ShowResult()
{
    auto tree = t.Tree;

    canvas c("PhotonCheck");
    for (const auto filter: {"","hitsCB == multiplicity",
                                "hitsTAPS == multiplicity",
                                "hitsTAPS > 0 && hitsCB > 0"})
    {
        c  << TTree_drawable(tree,"openings * 180 / 3.1415",filter)
           << drawoption("colz")
           << TTree_drawable(tree,"Theta() * 180 / 3.1415",filter) << endrow();
    }

    auto hits  = HistFac.makeTH1D("hits","","#",BinSettings(3),"hits");

    for (int i = 0 ; i < tree->GetEntries() ; ++i)
    {
        tree->GetEntry(i);
        if(t.hitsCB == t.multiplicity)
            hits->Fill("CB",1);
        if(t.hitsTAPS == t.multiplicity)
            hits->Fill("TAPS",1);
        if(t.hitsCB > 0 && t.hitsTAPS > 0 )
            hits->Fill("CB and TAPS",1);
    }

  c
            << hits
            << drawoption("colz")
            << TTree_drawable(tree,"hitsTAPS:hitsCB")
     << endc;

}




AUTO_REGISTER_PHYSICS(MCPhotonPairCheck)
