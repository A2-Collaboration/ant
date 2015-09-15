#pragma once

#include "analysis/physics/Physics.h"

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

    struct list_of_hists_t {
        virtual std::list<TH1*> Hists() =0;
    };

    struct candidatesEvent_t: list_of_hists_t {
        TH1D* nPerEvent;
        TH2D* nPerEventPerE;
        TH1D* splitPerEvent;
        TH2D* splitPos;
        candidatesEvent_t(SmartHistFactory& f, const std::string& prefix);
        void Fill(const data::ParticlePtr& mctrue, const data::CandidateList& cand);

        std::list<TH1*> Hists() override;
        TH3D* multiplicity_map;

        std::vector<TH1*> mult2_split_angles;

        TH1D* nCharged;
        TH2D* cluserSize_TAPS;
        TH2D* cluserSize_CB;
        TH2D* dEE_TAPS;
        TH2D* dEE_CB;
    };

    candidatesEvent_t nPerEvent;
    std::map<const ParticleTypeDatabase::Type*, candidatesEvent_t> nPerEvent_type;

public:
    ReconstructCheck(PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
