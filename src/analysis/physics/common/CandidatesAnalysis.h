#pragma once

#include "analysis/physics/Physics.h"
#include "utils/ClusterTools.h"

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
    TH2D* psa = nullptr;
    TH2D* psa_all = nullptr;
    TH2D* psa_all_angles = nullptr;
    TH1D* lateral_moment_cb = nullptr;
    TH1D* lateral_moment_taps = nullptr;


    utils::ClusterTools clustertools;

public:
    CandidatesAnalysis(const std::string& name,OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
