#include "ParticleID.h"
#include "data/Event.h"
#include "data/Candidate.h"
#include "data/Particle.h"
#include "TCutG.h"

using namespace ant;
using namespace ant::analysis;

BasicParticleID::~BasicParticleID()
{

}

void BasicParticleID::Process(const CandidatePtr& cand, ParticleList &particles_out) const
{
    ParticleTypeDatabase::Type* type = nullptr;
    bool hadronic(false);
    bool charged(false);

    if(tof) {
        hadronic = tof->IsInside(cand->ClusterEnergy(), cand->Time());
    }

    charged = cand->VetoEnergy() > 0.0;

    if(!charged) {
        if(hadronic) {
//            type=
        }
    }

}


CBTAPSBasicParticleID::~CBTAPSBasicParticleID()
{

}

void CBTAPSBasicParticleID::Process(const CandidatePtr& cand, ParticleList& particles_out) const
{
    if(cand->Detector() & Detector_t::Any_t::CB()) {
        cb.Process(cand, particles_out);
    } else if(cand->Detector() & Detector_t::Any_t::TAPS()) {
        taps.Process(cand, particles_out);
    }
}

void ParticleID::Process(const CandidateList& cands, ParticleList& particles_out) const
{
    for(const auto& cand : cands) {
        Process(cand, particles_out);
    }
}
