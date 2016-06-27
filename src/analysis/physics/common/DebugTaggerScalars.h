#pragma once

#include "analysis/physics/Physics.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


class DebugTaggerScalars: public Physics {
protected:
    unsigned seenEvents = 0;
    TH1D* byChannel;
public:
    DebugTaggerScalars(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~DebugTaggerScalars();

    virtual void ProcessEvent(const TEvent&, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
