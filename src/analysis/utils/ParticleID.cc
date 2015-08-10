#include "ParticleID.h"
#include "data/Event.h"
#include "data/Candidate.h"
#include "data/Particle.h"
#include "TCutG.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::utils;

BasicParticleID::~BasicParticleID()
{

}

bool TestCut(const std::shared_ptr<TCutG>& cut, const double& x, const double& y) {
    return (cut) && cut->IsInside(x,y);
}

const ParticleTypeDatabase::Type* BasicParticleID::Identify(const std::shared_ptr<ant::Candidate>& cand) const
{
    const bool hadronic =    TestCut(tof,  cand->ClusterEnergy(), cand->Time())
                  || TestCut(size, cand->ClusterEnergy(), cand->ClusterSize());

    const bool hadronic_enabled = (tof) || (size);

    const bool charged = cand->VetoEnergy() > 0.0;



    if(!charged) {

        if(hadronic) {
            return addressof(ParticleTypeDatabase::Neutron);
        } else {
            return addressof(ParticleTypeDatabase::Photon);
        }

    } else {  //charged

        if(
           (hadronic_enabled && hadronic)
           || (TestCut(dEE_proton, cand->ClusterEnergy(), cand->VetoEnergy()))
           ) {
            return addressof(ParticleTypeDatabase::Proton);
        }

        if(
           TestCut(dEE_pion, cand->ClusterEnergy(), cand->VetoEnergy())
           ) {
            return addressof(ParticleTypeDatabase::PiCharged);
        }

        if(
           TestCut(dEE_electron, cand->ClusterEnergy(), cand->VetoEnergy())
           ) {
            return addressof(ParticleTypeDatabase::eCharged);
        }
    }

    return nullptr;
}

std::shared_ptr<Particle> BasicParticleID::Process(std::shared_ptr<ant::Candidate>& cand) const
{
    auto type = Identify(cand);
    if(type !=nullptr) {
       return std::make_shared<Particle>(*type, cand);
    }

    return nullptr;
}



CBTAPSBasicParticleID::~CBTAPSBasicParticleID()
{

}

std::shared_ptr<Particle> CBTAPSBasicParticleID::Process(std::shared_ptr<Candidate>& cand) const
{
    if(cand->Detector() & Detector_t::Any_t::CB) {
        return cb.Process(cand);
    } else if(cand->Detector() & Detector_t::Any_t::TAPS) {
        return taps.Process(cand);
    }

    return nullptr;
}
