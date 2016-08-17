#pragma once

#include "analysis/physics/Physics.h"
#include "utils/TimeSmearingHack.h"

#include <string>

class TH1D;
class TH2D;
class TH3D;


namespace ant {

struct hstack;

namespace analysis {
namespace physics {

class MCReconstructCheck : public Physics {
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
        PositionMapCB(HistogramFactory& f, const std::string& name, const std::string &title="");
        virtual void Fill(const double theta, const double phi, const double v=1.0) override;
        virtual ~PositionMapCB() = default;
    };

    struct PositionMapTAPS : PositionMap {
        PositionMapTAPS(HistogramFactory& f, const std::string& name, const std::string &title="");
        virtual void Fill(const double cos, const double phi, const double v=1.0) override;
        virtual void Draw(const std::string& option) const override;
        virtual ~PositionMapTAPS() = default;
    };



    struct histgroup {
        const std::string Prefix;

        TH1D* nPerEvent;
        TH2D* nPerEventPerE;
        TH1D* splitPerEvent;
        std::shared_ptr<PositionMap> splitFlagPos;
        std::shared_ptr<PositionMap> splitPos;
        std::shared_ptr<PositionMap> touchesholeFlagPos;

        std::vector<TH1D*> mult2_split_angles;
        ant::hstack* splitstack;

        TH1D* nCharged;
        TH2D* cluserSize;
        TH2D* cluserSize_true;
        TH2D* dEE;
        TH2D* dEE_true;
        std::shared_ptr<PositionMap> posCharged;
        TH1D* unmatched_veto;
        TH1D* veto_cand_phi_diff;

        TH2D* energyinout;

        TH2D* thetainout;
        TH2D* phiinout;
        TH2D* anglediff;

        std::shared_ptr<PositionMap> mult1_positions;
        std::shared_ptr<PositionMap> energy_recov;

        std::shared_ptr<PositionMap> input_positions;
        std::shared_ptr<PositionMap> mult1_chargedPos;



        enum class detectortype {
            All, CB, TAPS
        };

        std::shared_ptr<PositionMap> makePosMap(HistogramFactory& f, detectortype d, const std::string& name, const std::string title="");

        histgroup(const HistogramFactory& parent, const std::string& prefix, detectortype d=detectortype::All);
        void Fill(const TParticlePtr& mctrue, const TCandidateList& cand, const TClusterList& all_clusters);
        void ShowResult() const;
        void Finish();

        histgroup(const histgroup&) = delete;
        histgroup& operator =(const histgroup&) = delete;

    };

    struct TAPSVetoMatch {
        TH2D* vetoElement_dist;
        TAPSVetoMatch(HistogramFactory& f);
        TAPSVetoMatch(const TAPSVetoMatch&) = delete;
        TAPSVetoMatch& operator =(const TAPSVetoMatch&) = delete;

        void ShowResult();
        void Finish() {}
        void Fill(const TCandidateList& cands, const TClusterList& all_clusters);
    };

    histgroup cb_group;
    histgroup taps_group;
    histgroup all_group;
    TAPSVetoMatch tapsveto;

    utils::TimeSmearingHack timesmear;


    const bool mult1_only;

    TTree* tree = nullptr;

    unsigned b_mult = 0;

    double b_rE     = 0.0;
    double b_rTheta = 0.0;
    double b_rPhi   = 0.0;
    double b_rVeto  = 0.0;
    double b_rTime  = 0.0;
    unsigned b_rSize = 0;

    double b_tE = 0.0;
    double b_tTheta = 0.0;
    double b_tPhi = 0.0;
    unsigned b_Cal = 0;

public:
    MCReconstructCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
