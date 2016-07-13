#include "analysis/physics/Physics.h"


namespace ant {
namespace analysis {
namespace physics {

class InterpolatedPulls : public Physics {


public:
    InterpolatedPulls(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;
};


}}}
