#ifndef ANTPHYSICS_H
#define ANTPHYSICS_H

#include "Event.h"
#include "plot/HistogramFactories.h"
#include <list>
#include <string>
#include <memory>

class TFile;
class TDirectory;

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

class PhysicsManager {
protected:
    using physics_list_t = std::list< std::unique_ptr<ant::Physics> >;

    physics_list_t physics;

public:
    PhysicsManager();
    template <typename T, typename ... args_t>
    void AddPhysics(args_t&&... args) {
        physics.emplace_back(
                    std::unique_ptr<ant::Physics>(
                        new T(std::forward<args_t>(args)...)));
    }
};
}


#endif
