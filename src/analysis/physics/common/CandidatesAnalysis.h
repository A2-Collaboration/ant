#pragma once

#include "analysis/physics/Physics.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class CandidatesAnalysis : public Physics {
protected:
    TH1D* nCandidatesEvent = nullptr;
    TH1D* CandMultiplicities = nullptr;
    TH1D* energy = nullptr;
    TH1D* theta = nullptr;
    TH1D* phi = nullptr;
    TH1D* ggIM = nullptr;
    TH1D* ttIM = nullptr;
    TH2D* cbdEE = nullptr;
    TH2D* cbtof = nullptr;
    TH2D* tapsdEE = nullptr;
    TH2D* tapstof = nullptr;
    TH1D* detectors = nullptr;

public:
    CandidatesAnalysis(const std::string& name,PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
