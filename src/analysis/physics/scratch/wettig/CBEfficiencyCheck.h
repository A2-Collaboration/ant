#pragma once

#include "analysis/physics/Physics.h"

#include <string>
#include <memory>
#include <vector>

class TH1D;

namespace ant {

class TH2CB;

namespace analysis {
namespace physics {


class CBEfficiencyCheck : public Physics {

protected:

    std::unique_ptr<TH2CB> hist_hits = nullptr;
    std::unique_ptr<TH2CB> hist_centClust = nullptr;
    std::unique_ptr<TH2CB> grid = nullptr;
    TCanvas* canv = nullptr;

    const int w_px;


public:
    CBEfficiencyCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    virtual ~CBEfficiencyCheck(){}
};

}
}
}
