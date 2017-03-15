#include "Omega_EpEm_mc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Omega_EpEm_mc::Omega_EpEm_mc(const string &name, OptionsPtr opts) :
    Physics(name, opts)
{
    BinSettings bins_nParticles(10);
    BinSettings bins_IM(200);

    // some tree
    t.CreateBranches(HistFac.makeTTree("t"));

    // histograms
    h_IM_2e = HistFac.makeTH1D("e^+ e^- IM",
                                  "M_{e^+ e^-} [MeV]","#",
                                  bins_IM,
                                  "IM_2e"
                                  );
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

    // loop over eCharged
    for (const auto& gE: eCharged)
    {
        t.p().emplace_back(*gE);
        const auto detHit = a2geo.DetectorFromAngles(gE->Theta(),
                                                     gE->Phi()    );
        if (detHit == Detector_t::Type_t::CB)
            t.hitsCB++;
        if (detHit == Detector_t::Type_t::TAPS)
            t.hitsTAPS++;
    }
    t.nEcharged = eCharged.size();

t.Tree->Fill(); // DO NOT FORGET
}

void Omega_EpEm_mc::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
        << h_IM_2e
        << TTree_drawable(t.Tree,"nEcharged")
        << TTree_drawable(t.Tree,"nEplus")
        << TTree_drawable(t.Tree,"nEminus")
        //<< TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
        << endc; // actually draws the canvas
}

void Omega_EpEm_mc::Finish()
{

}


AUTO_REGISTER_PHYSICS(Omega_EpEm_mc)
