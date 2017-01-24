#pragma once

#include "physics/Physics.h"

#include <map>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class TestParticleCombinatorics: public Physics {
protected:
    TH1D* ggim;
    TH1D* gggim;
    TH1D* nphotons;
    TH1D* nprotons;

    std::map<const ant::ParticleTypeDatabase::Type*, TH1D*> EHists;

public:
    TestParticleCombinatorics(const std::string& name, OptionsPtr opts);
    virtual ~TestParticleCombinatorics() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
