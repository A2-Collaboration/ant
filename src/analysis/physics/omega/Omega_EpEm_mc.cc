#include "Omega_EpEm_mc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Omega_EpEm_mc::Omega_EpEm_mc(const string &name, OptionsPtr opts) :
    Physics(name, opts)
{
    // some tree
    t.CreateBranches(HistFac.makeTTree("t"));
}

void Omega_EpEm_mc::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& particleTree   = event.MCTrue().ParticleTree;
    if(!particleTree)
        return;

    // get electrons/positrons;
    vector<TParticlePtr> eCharged;
    particleTree->Map_nodes([&eCharged](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::eCharged)
            eCharged.emplace_back(pt->Get());
    });

}

void Omega_EpEm_mc::ShowResult()
{

}

void Omega_EpEm_mc::Finish()
{

}


AUTO_REGISTER_PHYSICS(Omega_EpEm_mc)
