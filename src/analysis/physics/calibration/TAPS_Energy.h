#pragma once

#include "analysis/physics/Physics.h"

#include "root-addons/cbtaps_display/TH2TAPS.h"


namespace ant {

namespace expconfig {
namespace detector {
struct TAPS;
}}

namespace analysis {
namespace physics {

class TAPS_Energy : public Physics {

protected:
    TH2D* ggIM_CBTAPS_PIDVetos = nullptr;
    TH2D* ggIM_CBTAPS_PID = nullptr;
    TH2D* ggIM_Any_PID = nullptr;
    TH2D* ggIM_Any_PIDVetos = nullptr;
    TH3D* ggIM_mult = nullptr;
    TH2D* timing_cuts = nullptr;
    TH2D* h_pedestals = nullptr;
    TH2TAPS* h_tapsdisplay = nullptr;

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
        void Fill(const TCandidate& cand);
    };

    TTree* cands_tree = nullptr;
    tree_data_t cands_CB;
    tree_data_t cands_TAPS;

public:

    TAPS_Energy(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics
