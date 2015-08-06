#pragma once

#include "data/Event.h"

#include <memory>

class TCutG;

namespace ant {
namespace analysis {

class ParticleID {
public:
    virtual ~ParticleID() {}

    /**
     * @brief Identify particles in this event.
     * @param event Is modified, particles added
     */
    virtual void Process(const ant::CandidateList& cands, ant::ParticleList& particles_out) const;

    virtual void Process(const CandidatePtr& cand, ParticleList &particles_out) const =0;
};


class BasicParticleID: public ParticleID {
public:
    virtual ~BasicParticleID();

    std::shared_ptr<TCutG> dEE_proton;
    std::shared_ptr<TCutG> dEE_pion;
    std::shared_ptr<TCutG> dEE_electron;

    std::shared_ptr<TCutG> tof;

    std::shared_ptr<TCutG> size;

    virtual void Process(const CandidatePtr& cand, ParticleList &particles_out) const override;
};

class CBTAPSBasicParticleID: public ParticleID {
protected:
    BasicParticleID cb;
    BasicParticleID taps;

public:
    virtual ~CBTAPSBasicParticleID();

    virtual void Process(const CandidatePtr& cand, ParticleList &particles_out) const override;
};

}
}

