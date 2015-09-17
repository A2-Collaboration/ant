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

    struct PositionMap: ant::root_drawable_traits {
        TH2D* maphist = nullptr;
        virtual void Fill(const double theta, const double phi, const double v=1.0) =0;
        virtual void Draw(const std::string& option) const override;

        virtual ~PositionMap() = default;

        PositionMap() = default;
        PositionMap(const PositionMap&) = delete;
        PositionMap& operator =(const PositionMap&) = delete;
    };

    struct PositionMapCB : PositionMap {
        PositionMapCB(SmartHistFactory& f, const std::string& name, const std::string &title="");
        virtual void Fill(const double theta, const double phi, const double v=1.0) override;
        virtual ~PositionMapCB() = default;
    };

    struct PositionMapTAPS : PositionMap {
        PositionMapTAPS(SmartHistFactory& f, const std::string& name, const std::string &title="");
        virtual void Fill(const double cos, const double phi, const double v=1.0) override;
        virtual void Draw(const std::string& option) const override;
        virtual ~PositionMapTAPS() = default;
    };



    struct histgroup {
        const std::string Prefix;

        TH1D* nPerEvent;
        TH2D* nPerEventPerE;
        TH1D* splitPerEvent;
        std::unique_ptr<PositionMap> splitPos;

        TH3D* multiplicity_map;

        std::vector<TH1D*> mult2_split_angles;

        TH1D* nCharged;
        TH2D* cluserSize;
        TH2D* dEE;
        TH2D* dEE_true;
        std::unique_ptr<PositionMap> posCharged;
        TH1D* unmatched_veto;
        TH2D* edge_flag_pos;

        TH1D* veto_cand_phi_diff;

        std::unique_ptr<PositionMap> energy_recov_norm;
        std::unique_ptr<PositionMap> energy_recov;

        enum class detectortype {
            All, CB, TAPS
        };

        std::unique_ptr<PositionMap> makePosMap(SmartHistFactory& f, detectortype d, const std::string& name, const std::string title="");

        histgroup(SmartHistFactory& f, const std::string& prefix, detectortype d=detectortype::All);
        void Fill(const data::ParticlePtr& mctrue, const data::CandidateList& cand, const data::ClusterList& insane);
        void ShowResult() const;
        void Finish();

        histgroup(const histgroup&) = delete;
        histgroup& operator =(const histgroup&) = delete;

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
