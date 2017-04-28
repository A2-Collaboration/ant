#pragma once

#include "analysis/physics/Physics.h"
#include "root-addons/analysis_codes/hstack.h"

namespace ant {
namespace analysis {
namespace physics {



class MCClusteringCheck : public Physics {
protected:

    const Detector_t::Type_t detectorType;

    struct opening_angle_t {

        const interval<double> opening_angle_range;
        opening_angle_t(const interval<double> opening_angle_range_,
                        const HistogramFactory& HistFac,
                        Detector_t::Type_t detectorType);

        TH1D* h_nCands = nullptr;
        TH1D* h_nSplits = nullptr;
        TH2D* h_ErecEtrue1 = nullptr;
        TH2D* h_ErecEtrue2 = nullptr;
        TH2D* h_OpeningAngle1 = nullptr;
        TH2D* h_OpeningAngle2 = nullptr;
        TH2D* h_DiffAngleTheta1 = nullptr;
        TH2D* h_DiffAngleTheta2 = nullptr;
        TH2D* h_DiffAnglePhi1 = nullptr;
        TH2D* h_DiffAnglePhi2 = nullptr;
        TH1D* h_nUnmatchedCandsMinAngle = nullptr;

        bool Fill(double opening_angle, unsigned nCands,
                  const TParticle& true_photon1, const TParticle& true_photon2,
                  const TCandidatePtr& best_cand1, const TCandidatePtr& best_cand2,
                  const TCandidatePtrList& unmatched_cands
                  ) const;
        void Show(ant::canvas& c) const;

    };

    const std::vector<opening_angle_t> opening_angles;

    TH1D* h_Steps = nullptr;
    std::vector<TH1D*> h_Cands_OpAng;

public:
    MCClusteringCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
