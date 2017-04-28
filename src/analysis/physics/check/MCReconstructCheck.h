#pragma once

#include "analysis/physics/Physics.h"
#include "utils/TimeSmearingHack.h"
#include "base/WrapTTree.h"

#include <string>

class TH1D;
class TH2D;
class TH3D;


namespace ant {

struct hstack;

namespace expconfig {
namespace detector {
struct TAPS;
}
}

namespace analysis {
namespace physics {

class MCReconstructCheck : public Physics {
protected:

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
        PositionMapCB(const HistogramFactory& f, const std::string& name, const std::string &title="");
        virtual void Fill(const double theta, const double phi, const double v=1.0) override;
        virtual ~PositionMapCB() = default;
    };

    struct PositionMapTAPS : PositionMap {
        const std::shared_ptr<const expconfig::detector::TAPS> taps;
        const double taps_dist;
        PositionMapTAPS(const HistogramFactory& f, const std::string& name, const std::string &title="");
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
        ant::hstack* mult2_split_stack;

        TH1D* nCharged;
        TH2D* cluserSize;
        TH2D* cluserSize_true;
        TH2D* dEE;
        TH2D* dEE_true;
        std::shared_ptr<PositionMap> posCharged;
        TH1D* unmatched_veto;
        TH1D* veto_cand_phi_diff;

        TH2D* energy_rec_true;

        TH2D* theta_diff_truetheta;
        TH2D* phi_diff_truetheta;
        TH2D* openingangle_diff_truetheta;

        TH2D* theta_diff_trueEk;
        TH2D* phi_diff_trueEk;
        TH2D* openingangle_diff_trueEk;

        std::shared_ptr<PositionMap> rec_mult1_positions;
        std::shared_ptr<PositionMap> mc_mult1_positions;
        std::shared_ptr<PositionMap> energy_recov;

        std::shared_ptr<PositionMap> input_positions;
        std::shared_ptr<PositionMap> mult1_chargedPos;



        enum class detectortype {
            All, CB, TAPS
        };

        std::shared_ptr<PositionMap> makePosMap(const HistogramFactory& f, detectortype d, const std::string& name, const std::string title="");

        histgroup(const HistogramFactory& parent, const std::string& prefix, detectortype d=detectortype::All);
        void Fill(const TParticlePtr& mctrue, const TCandidateList& cand, const TClusterList& all_clusters);
        void ShowResult() const;
        void Finish();

        histgroup(const histgroup&) = delete;
        histgroup& operator =(const histgroup&) = delete;
    };

    struct TAPSVetoMatch {
        TH2D* vetoElement_dist;
        TAPSVetoMatch(const HistogramFactory& f);
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
    const bool show_cb_only;
    const bool fill_tree;


    struct tree_t : WrapTTree {
        ADD_BRANCH_T(unsigned, mult)

        // reconstructed variables
        ADD_BRANCH_T(double,   rE )
        ADD_BRANCH_T(double,   rTheta)
        ADD_BRANCH_T(double,   rPhi)
        ADD_BRANCH_T(double,   rVeto)
        ADD_BRANCH_T(double,   rTime)
        ADD_BRANCH_T(unsigned, rSize)

        // true variables
        ADD_BRANCH_T(double,   tE)
        ADD_BRANCH_T(double,   tTheta)
        ADD_BRANCH_T(double,   tPhi)

        // Cal=1 is CB, Cal=2 is TAPS
        ADD_BRANCH_T(unsigned, Cal)
    };

    tree_t t;


public:
    MCReconstructCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
