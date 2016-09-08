#include "analysis/physics/Physics.h"

#include "analysis/plot/PromptRandomHist.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/FitterSergey.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class EtapSergey : public Physics {
public:

protected:

    utils::FitterSergey fitter_sergey;

public:
    EtapSergey(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}
