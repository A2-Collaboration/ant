#ifndef ANTPHYSICS_H
#define ANTPHYSICS_H

#include "Event.h"
#include "plot/HistogramFactories.h"

namespace ant {

class Physics {
protected:
    SmartHistFactory HistFac;
public:
    Physics(const std::string& name);
    virtual ~Physics() {}
    virtual void ProcessEvent(const ant::Event& event) =0;
    virtual void Finish() =0;
    virtual void ShowResult() =0;
};

class DebugPhysics: public Physics {
public:
    DebugPhysics(const std::string& name="DebugPhysics"): Physics(name) {}
    virtual ~DebugPhysics() {}

    virtual void ProcessEvent(const ant::Event& event);
    virtual void Finish();
    virtual void ShowResult();
};
}


#endif
