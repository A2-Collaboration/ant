#include "Omega_EpEm_mc.h"

//#include "base/Logger.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Omega_EpEm_mc::Omega_EpEm_mc(const string &name, OptionsPtr opts) :
    Physics(name, opts)
{
//    BinSettings bins_nParticles(10);
//    BinSettings bins_IM(200);

    // create the tree which contains all information
    t.CreateBranches(HistFac.makeTTree("t"));

}


void Omega_EpEm_mc::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& particleTree   = event.MCTrue().ParticleTree;
    if(!particleTree)
        return;

    // get electrons/positrons;
    vector<TParticlePtr> eCharged;
    vector<TParticlePtr> ePlus;
    vector<TParticlePtr> eMinus;
    vector<TParticlePtr> proton;

    particleTree->Map_nodes([&eCharged](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::eCharged)
            eCharged.emplace_back(pt->Get());
    });
    particleTree->Map_nodes([&ePlus](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::ePlus)
            ePlus.emplace_back(pt->Get());
    });
    particleTree->Map_nodes([&eMinus](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::eMinus)
            eMinus.emplace_back(pt->Get());
    });
    particleTree->Map_nodes([&proton](const TParticleTree_t& pt){
        if (pt->Get()->Type() == ParticleTypeDatabase::Proton)
            proton.emplace_back(pt->Get());
    });

    // loop over eCharged
    for (const auto& gE: eCharged)
    {
        t.eVector().emplace_back(*gE);
        const auto detHit = a2geo.DetectorFromAngles(gE->Theta(),
                                                     gE->Phi()    );
        if (detHit == Detector_t::Type_t::CB)
            t.hitsCB++;
        if (detHit == Detector_t::Type_t::TAPS)
            t.hitsTAPS++;
    }

    // loop over proton (should be 1 anyway)
    for (const auto& gE: proton)
    {
        const auto detHit = a2geo.DetectorFromAngles(gE->Theta(),
                                                     gE->Phi()    );
        if (detHit == Detector_t::Type_t::CB)
            t.hitsCB++;
        if (detHit == Detector_t::Type_t::TAPS)
            t.hitsTAPS++;
    }


    if (eCharged.size() == 2)
    {
        auto combs = utils::makeCombination(eCharged,2);
        t.eeOpenAngle = TParticle::CalcAngle(combs.at(0), combs.at(1));
        const auto& c1 = combs.at(0);
        const auto& c2 = combs.at(1);
        const auto sum = (*c1 + *c2);
        t.eeIM = sum.M();
        const auto& boost_e1 = Boost(*c1, -sum.BoostVector());
        const auto& boost_e2 = Boost(*c2, -sum.BoostVector());
        t.eeBoostOpenAngle = boost_e1.Angle(boost_e2);
    }

    if (proton.size() == 1)
    {
        t.pEk = proton.at(0)->Ek();
        t.pTheta = proton.at(0)->Theta();
        t.pPhi = proton.at(0)->Phi();
    }
    if (eMinus.size() == 1)
    {
        t.eMEk = eMinus.at(0)->Ek();
        t.eMTheta = eMinus.at(0)->Theta();
        t.eMPhi = eMinus.at(0)->Phi();
    }
    if (ePlus.size() == 1)
    {
        t.ePEk = ePlus.at(0)->Ek();
        t.ePTheta = ePlus.at(0)->Theta();
        t.ePPhi = ePlus.at(0)->Phi();
    }

t.fillAndReset(); // do not forget!
}

void Omega_EpEm_mc::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
//        << drawoption("colz") << TTree_drawable(t.Tree,"hitsCB:hitsTAPS")
//        << TTree_drawable(t.Tree,"hitsCB")
//        << TTree_drawable(t.Tree,"hitsTAPS")
        << TTree_drawable(t.Tree,"eeOpenAngle * 180 / 3.1415")
        << TTree_drawable(t.Tree,"eeBoostOpenAngle * 180 / 3.1415")
        << TTree_drawable(t.Tree,"eeIM")
        << drawoption("colz") << TTree_drawable(t.Tree,"eMTheta * 180 / 3.1415:eMEk")
        << drawoption("colz") << TTree_drawable(t.Tree,"ePTheta * 180 / 3.1415:ePEk")
        << drawoption("colz") << TTree_drawable(t.Tree,"pTheta * 180 / 3.1415:pEk")
        << TTree_drawable(t.Tree,"eMPhi")
        << TTree_drawable(t.Tree,"ePPhi")
        << TTree_drawable(t.Tree,"pPhi")
        << endc; // actually draws the canvas
}

void Omega_EpEm_mc::Finish()
{

}


AUTO_REGISTER_PHYSICS(Omega_EpEm_mc)
