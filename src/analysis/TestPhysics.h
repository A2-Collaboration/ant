#ifndef TESTPHYSICS_H
#define TESTPHYSICS_H

#include "AntPhysics.h"
#include "plot/Histogram.h"
#include "plot/SmartHist.h"

#include <map>

namespace ant {
class ParticleCombinatoricsTest: public ant::Physics {
protected:
    SmartHist1<double> ggim;
    SmartHist1<double> gggim;
    SmartHist1<int>    nphotons;
    SmartHist1<int>    nprotons;

    std::map<const ant::ParticleTypeDatabase::Type*, SmartHist1<const ParticlePtr&>> EHists;


public:
    ParticleCombinatoricsTest(const std::string& name="ParticleCombinatoricsTest");
    virtual ~ParticleCombinatoricsTest() {}

    virtual void ProcessEvent(const ant::Event& event);
    virtual void Finish();
    virtual void ShowResult();
};

}
#endif
