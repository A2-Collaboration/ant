#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"
#include "plot/PromptRandomHist.h"
#include "utils/KinFitter.h"
#include "expconfig/detectors/TAPS.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {
namespace physics {

class EtapProton : public Physics {

protected:
    PiecewiseInterval<unsigned> multiplicities;

    TTree* tree = nullptr;

    unsigned  b_nCB = 0;
    unsigned  b_nTAPS = 0;
    double    b_CBAvgTime  = 0.0;
    double    b_CBSumVetoE  = 0.0;

    TCandidate b_Proton;
    std::vector<TCandidate> b_Photons;
    TLorentzVector b_PhotonSum;
    double b_ProtonCopl;
    double b_ProtonBeta;

    PromptRandom::Switch promptrandom;

    std::vector<std::unique_ptr<utils::KinFitter>> fitters;
    double b_BestChi2;
    unsigned b_NGoodFits;
    unsigned b_NFitIterations;
    double b_TaggW;
    double b_TaggE;
    double b_TaggT;
    unsigned b_TaggCh;
    unsigned b_TaggN;

    TLorentzVector b_FittedProton;
    std::vector<TLorentzVector> b_FittedPhotons;
    TLorentzVector b_FittedPhotonSum;
    double b_FittedProtonCopl;

    std::shared_ptr<expconfig::detector::TAPS> taps_detector;

public:

    EtapProton(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
