#pragma once

#include "physics/Physics.h"
#include "plot/Histogram.h"
#include "plot/SmartHist.h"

#include <map>

namespace ant {
namespace analysis {
namespace physics {

class ParticleCombinatoricsTest: public Physics {
protected:
    SmartHist1<double> ggim;
    SmartHist1<double> gggim;
    SmartHist1<int>    nphotons;
    SmartHist1<int>    nprotons;

    std::map<const ant::ParticleTypeDatabase::Type*, SmartHist1<const data::ParticlePtr&>> EHists;

public:
    ParticleCombinatoricsTest(const std::string& name, PhysOptPtr opts);
    virtual ~ParticleCombinatoricsTest() {}

    virtual void ProcessEvent(const data::Event& event);
    virtual void Finish();
    virtual void ShowResult();
};

}
}
}
