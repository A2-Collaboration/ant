#pragma once

#include "analysis/physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include <vector>

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class IMPlots : public Physics {
public:

    PromptRandom::Switch prs;
    std::vector<PromptRandom::Hist1> m;
    unsigned MinNGamma() const noexcept { return 2;}
    unsigned MaxNGamma() const noexcept { return unsigned(m.size())+2; }




public:
    IMPlots(const std::string& name, PhysOptPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


class Symmetric2Gamma : public Physics {
protected:
    TH1D* h_symIM = nullptr;
    TTree* tree   = nullptr;

    double b_IM;
    double b_E;
    double b_E1;
    double b_E2;
    double b_theta1;
    double b_theta2;
    double b_phi1;
    double b_phi2;

    double perc = 0.2;

public:
    Symmetric2Gamma(const std::string& name, PhysOptPtr opts);
    virtual ~Symmetric2Gamma();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
