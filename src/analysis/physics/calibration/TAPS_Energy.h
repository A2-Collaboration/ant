#pragma once

#include "analysis/physics/Physics.h"

namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace analysis {
namespace physics {

class TAPS_Energy : public Physics {

protected:
    TH2D* ggIM = nullptr;
    TH3D* ggIM_mult = nullptr;
    TH2D* timing_cuts = nullptr;
    TH2D* h_pedestals = nullptr;

    std::shared_ptr<expconfig::detector::TAPS> taps_detector;

    struct tree_data_t {
        std::vector<double> Ek;
        std::vector<double> Theta;
        std::vector<double> Phi;
        std::vector<double> VetoE;
        std::vector<double> Time;
        std::vector<unsigned> Channel;
        void Setup(const std::string& prefix, TTree* tree);
        void Clear();
        void Fill(const data::Candidate& cand);
    };

    TTree* cands_tree = nullptr;
    tree_data_t cands_CB;
    tree_data_t cands_TAPS;

public:

    TAPS_Energy(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
};

}}} // namespace ant::analysis::physics