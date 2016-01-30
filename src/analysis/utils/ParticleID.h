#pragma once

#include "tree/TParticle.h"

#include "base/ParticleType.h"

#include <memory>

class TCutG;


namespace ant {
class WrapTFile;
namespace analysis {


namespace utils {

class ParticleID {
public:
    virtual ~ParticleID() {}

    virtual const ParticleTypeDatabase::Type* Identify(const TCandidatePtr& cand) const =0;

    virtual TParticlePtr Process(const TCandidatePtr& cand) const;
};


class SimpleParticleID: public ParticleID {
public:
    SimpleParticleID() {}
    virtual ~SimpleParticleID() {}

    virtual const ParticleTypeDatabase::Type* Identify(const TCandidatePtr& cand) const override;
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

    virtual const ParticleTypeDatabase::Type* Identify(const TCandidatePtr& cand) const override;
};

class CBTAPSBasicParticleID: public ParticleID {
protected:
    BasicParticleID cb;
    BasicParticleID taps;
    virtual void LoadFrom(WrapTFile& file);

public:
    CBTAPSBasicParticleID(const std::string& pidcutsdir);
    virtual ~CBTAPSBasicParticleID();

    virtual const ParticleTypeDatabase::Type* Identify(const TCandidatePtr& cand) const override;
};

}
}
}
