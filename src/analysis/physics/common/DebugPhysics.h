#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {


class DebugPhysics: public Physics {
protected:
    const bool noDump;
    const unsigned writeEvents;
    const bool keepReadHits;
    const bool requestSlowControl;
    unsigned seenEvents = 0;
    TID lastTID;
public:
    DebugPhysics(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~DebugPhysics();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

class DebugPIDAlignment: public Physics {
protected:
    TH2D* angles_mc = nullptr;
    TH2D* angles_candidates = nullptr;
    TH2D* angles_clusters = nullptr;
    TH1D* angles_diff = nullptr;
    TH1D* angles_diff_wrap = nullptr;
public:
    DebugPIDAlignment(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~DebugPIDAlignment();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
