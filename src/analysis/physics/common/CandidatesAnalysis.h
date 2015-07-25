#pragma once

#include "analysis/physics/Physics.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {

class CandidatesAnalysis : public ant::Physics {
protected:
    TH1D* nCandidatesEvent = nullptr;
    TH1D* energy = nullptr;
    TH1D* theta = nullptr;
    TH1D* phi = nullptr;
    TH1D* ggIM = nullptr;
    TH1D* ttIM = nullptr;

public:
    CandidatesAnalysis(const std::string& name="CandidatesAnalysis");

    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
