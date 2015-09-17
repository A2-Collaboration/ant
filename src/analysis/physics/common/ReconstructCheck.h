#pragma once

#include "analysis/physics/Physics.h"
#include "data/Candidate.h"

#include <string>

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class ReconstructCheck : public Physics {
protected:

    TH2D* EnergyRec_cb;
    TH2D* EnergyRec_taps;

    struct histgroup {
        const std::string Prefix;

        TH1D* nPerEvent;
        TH2D* nPerEventPerE;
        TH1D* splitPerEvent;
        TH2D* splitPos;

        TH3D* multiplicity_map;

        std::vector<TH1D*> mult2_split_angles;

        TH1D* nCharged;
        TH2D* cluserSize;
        TH2D* dEE;
        TH2D* dEE_true;
        TH2D* posCharged;
        TH1D* unmatched_veto;
        TH2D* edge_flag_pos;

        TH1D* veto_cand_phi_diff;

        TH2D* energy_recov_norm;
        TH2D* energy_recov;

        histgroup(SmartHistFactory& f, const std::string& prefix);
        void Fill(const data::ParticlePtr& mctrue, const data::CandidateList& cand, const data::ClusterList& insane);
        void ShowResult() const;
        void Finish();

    };

    histgroup cb_group;
    histgroup taps_group;
    histgroup all_group;

public:
    ReconstructCheck(PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
