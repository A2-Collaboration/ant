#pragma once

#include "ParticleType.h"

#include <memory>

class TCutG;


namespace ant {

class Candidate;
class Particle;

namespace analysis {

class ParticleID {
public:
    virtual ~ParticleID() {}

    virtual std::shared_ptr<Particle> Process(std::shared_ptr<ant::Candidate>& cand) const =0;
};


class BasicParticleID: public ParticleID {
public:
    BasicParticleID() {}
    virtual ~BasicParticleID();

    std::shared_ptr<TCutG> dEE_proton;
    std::shared_ptr<TCutG> dEE_pion;
    std::shared_ptr<TCutG> dEE_electron;

    std::shared_ptr<TCutG> tof;

    std::shared_ptr<TCutG> size;

    virtual const ParticleTypeDatabase::Type* Identify(const std::shared_ptr<Candidate>& cand) const;
    virtual std::shared_ptr<Particle> Process(std::shared_ptr<ant::Candidate>& cand) const override;
};

class CBTAPSBasicParticleID: public ParticleID {
protected:
    BasicParticleID cb;
    BasicParticleID taps;

public:
    CBTAPSBasicParticleID() {}
    virtual ~CBTAPSBasicParticleID();


    virtual std::shared_ptr<Particle> Process(std::shared_ptr<ant::Candidate>& cand) const override;
};

}
}

